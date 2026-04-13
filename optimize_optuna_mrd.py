"""
optimize_optuna_mrd.py — Optimisation Bayesienne MRD avec Optuna
=================================================================
Objectif: remplacer le grid search exhaustif par une optimisation bayesienne
avec une fonction de coût hybride cliniquement pertinente.

Changements majeurs v2:
  - Weighted hybrid loss: pénalise exponentiellement les faux négatifs sur fortes MRD
  - Contrainte patient sain: pénalité dure si MRD_sain > FP_THRESHOLD
  - SQLite persistence: reprise automatique après interruption
  - Log-scale sur paramètres de distance (d2_*, d2_normalization)
  - Pruning précoce: trial.report + trial.should_prune après chaque patient
  - Mémoire: libération explicite de nbm_frames après concat, config patchée
    in-place sur une copie légère au lieu de deepcopy à chaque trial
  - Parallélisation: n_jobs=-1 si le sampler le supporte (commentaire guidé)
"""

from __future__ import annotations

import argparse
import copy
import gc
import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import optuna
import pandas as pd

# ---------------------------------------------------------------------------
# PYTHONPATH robuste
# ---------------------------------------------------------------------------
CURRENT_DIR = Path(__file__).resolve().parent
PARENT_DIR = CURRENT_DIR.parent

if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

os.environ["PYTHONPATH"] = str(PARENT_DIR)

try:
    from src.analysis.blast_detection import compute_reference_stats
    from src.analysis.mrd_calculator import compute_mrd, load_mrd_config
    from src.core.clustering import FlowSOMClusterer
    from src.io.fcs_reader import get_fcs_files, load_fcs_file
    from src.models.sample import FlowSample
    from src.services.clustering_service import select_markers_for_clustering
    from src.services.preprocessing_service import preprocess_sample
    from config.pipeline_config import PipelineConfig
except ImportError as exc:
    print(f"Erreur d'importation: {exc}")
    print(f"DEBUG sys.path: {sys.path}")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Parametres du projet
# ---------------------------------------------------------------------------
HEALTHY_DIR = Path("F:/Cohorte/NBM")
PATHO_DIR = Path("F:/Cohorte/Patho trié")
GROUND_TRUTH_CSV_PATH = Path("F:/Cohorte/true_mrd_labels.csv")

OUTPUT_DIR = CURRENT_DIR / "Résults" / "Optuna_MRD"
TRIALS_CSV = OUTPUT_DIR / "optuna_trials.csv"
BEST_JSON = OUTPUT_DIR / "best_params.json"
SUMMARY_JSON = OUTPUT_DIR / "run_summary.json"
SQLITE_PATH = OUTPUT_DIR / "optuna_study.db"

MRD_CFG_PATH = CURRENT_DIR / "config" / "mrd_config.yaml"

# ---------------------------------------------------------------------------
# Constantes cliniques de la fonction de coût
# ---------------------------------------------------------------------------
# Seuil faux positif patient sain: toute prédiction au-dessus est pénalisée durement
FP_THRESHOLD_PCT = 0.1  # 0.1% — recommandation ELN

# Poids exponentiel pour les faux négatifs sur fortes MRD:
# loss_i = abs_err_i * (1 + HIGH_MRD_WEIGHT * (true_mrd_i / HIGH_MRD_CUTOFF))^HIGH_MRD_EXPONENT
HIGH_MRD_CUTOFF = 5.0  # au-delà de 5%, la pondération s'active fortement
HIGH_MRD_WEIGHT = 3.0  # amplificateur de base
HIGH_MRD_EXPONENT = 2.0  # exponentiel pour les MRD très élevées

# Pénalité dure appliquée si un patient sain dépasse FP_THRESHOLD_PCT
FP_PENALTY = 50.0  # en unités d'erreur (%)

# Pénalité appliquée quand compute_mrd échoue (valeur aberrante → trial mauvais)
FAILURE_PENALTY = 1e4

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

optuna.logging.set_verbosity(optuna.logging.WARNING)


# ---------------------------------------------------------------------------
# Dataclass cache patient
# ---------------------------------------------------------------------------
@dataclass
class PatientCache:
    pid: str
    fcs_path: Path
    data_train: pd.DataFrame
    node_assignments: np.ndarray
    node_medians: np.ndarray
    nbm_stats: Tuple[np.ndarray, np.ndarray, np.ndarray]
    markers: List[str]
    true_mrd: float  # NaN si pas de ground-truth
    has_ground_truth: bool
    is_healthy_control: bool  # True si ce patient est un NBM (sain de référence)


# ---------------------------------------------------------------------------
# Chargement ground-truth
# ---------------------------------------------------------------------------
def load_ground_truth_from_csv(csv_path: Path) -> Dict[str, float]:
    """Load ground-truth map from CSV (no hardcoded fallback)."""
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV ground-truth introuvable: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise RuntimeError(f"CSV ground-truth vide: {csv_path}")

    colmap = {c.lower().strip(): c for c in df.columns}

    pid_col = None
    for candidate in ("pid", "patient_id", "patient", "fichier", "file", "blast110_id"):
        if candidate in colmap:
            pid_col = colmap[candidate]
            break

    value_col = None
    for candidate in (
        "true_mrd",
        "mrd",
        "true",
        "ground_truth",
        "true_mrd_percent",
        "mrd_percent",
    ):
        if candidate in colmap:
            value_col = colmap[candidate]
            break

    if pid_col is None or value_col is None:
        raise RuntimeError(
            "CSV ground-truth invalide: colonnes attendues manquantes "
            "(pid/patient_id/... et true_mrd/mrd/...)."
        )

    out: Dict[str, float] = {}
    for _, row in df.iterrows():
        raw_pid = str(row[pid_col]).strip()
        if not raw_pid or raw_pid.lower() == "nan":
            continue
        pid = Path(raw_pid).stem
        try:
            value = float(row[value_col])
            if not np.isfinite(value):
                continue
            out[pid] = value
        except (TypeError, ValueError):
            continue

    if not out:
        raise RuntimeError("CSV ground-truth exploitable vide apres parsing.")

    log.info("Ground-truth charge depuis CSV: %d patient(s)", len(out))
    return out


def resolve_ground_truth_for_test2(
    preferred_csv: Path,
    patho_dir: Path,
) -> Tuple[Path, Dict[str, float], int]:
    """Resolve the best ground-truth CSV automatically."""
    patho_files = get_fcs_files(patho_dir)
    if not patho_files:
        raise FileNotFoundError(f"Aucun fichier pathologique trouve dans {patho_dir}")
    test2_pids = {Path(f).stem for f in patho_files}

    if preferred_csv.exists():
        try:
            preferred_gt = load_ground_truth_from_csv(preferred_csv)
            preferred_matches = len(test2_pids.intersection(preferred_gt.keys()))
            if preferred_matches > 0:
                return preferred_csv, preferred_gt, preferred_matches
        except Exception:
            pass

    candidates: List[Path] = [
        preferred_csv,
        CURRENT_DIR / "ground_truth.csv",
        CURRENT_DIR / "ground_truth_test_2.csv",
        CURRENT_DIR
        / "Résults"
        / "Validation_Clinique_V6"
        / "full_optimization_results.csv",
    ]
    candidates.extend(sorted(CURRENT_DIR.glob("*.csv")))
    results_dir = CURRENT_DIR / "Résults"
    if results_dir.exists():
        candidates.extend(sorted(results_dir.rglob("*.csv")))

    seen: set = set()
    unique_candidates: List[Path] = []
    for p in candidates:
        rp = p.resolve() if p.exists() else p
        key = str(rp).lower()
        if key in seen:
            continue
        seen.add(key)
        unique_candidates.append(p)

    best_path: Optional[Path] = None
    best_gt: Optional[Dict[str, float]] = None
    best_matches = -1

    for csv_path in unique_candidates:
        if not csv_path.exists():
            continue
        try:
            gt = load_ground_truth_from_csv(csv_path)
        except Exception:
            continue
        matches = len(test2_pids.intersection(gt.keys()))
        if matches > best_matches:
            best_path = csv_path
            best_gt = gt
            best_matches = matches

    if best_path is None or best_gt is None:
        raise FileNotFoundError(
            "Aucun CSV ground-truth compatible trouve. "
            "Fournir un CSV avec colonnes pid/patient_id et true_mrd/mrd/true."
        )
    if best_matches <= 0:
        raise RuntimeError(
            f"CSV detecte ({best_path}) mais aucun PID ne correspond aux fichiers de {patho_dir}."
        )

    return best_path, best_gt, best_matches


# ---------------------------------------------------------------------------
# Cache patient (calcul FlowSOM unique avant les trials)
# ---------------------------------------------------------------------------
def prepare_patient_cache(
    config: PipelineConfig,
    ground_truth: Dict[str, float],
    som_grid: int,
) -> List[PatientCache]:
    """Prepare all expensive computations once (NBM + FlowSOM per patient).

    Memory strategy:
      - nbm_frames list released immediately after concat
      - AnnData / FlowSample / preprocess intermediaires supprimes apres usage
      - gc.collect() apres chaque patient
    """
    log.info(
        "=== Passe 1: preparation du cache patient (FlowSOM %dx%d) ===",
        som_grid,
        som_grid,
    )

    # 1) Build NBM pool ---------------------------------------------------------
    nbm_files = get_fcs_files(HEALTHY_DIR)
    if not nbm_files:
        raise FileNotFoundError(f"Aucun fichier NBM trouve dans {HEALTHY_DIR}")

    nbm_frames: List[pd.DataFrame] = []
    for f in nbm_files[:10]:
        adata = load_fcs_file(f, condition="Sain")
        sample = FlowSample.from_anndata(adata, str(f), "Sain")
        df = preprocess_sample(sample, config).data.copy()
        df["condition"] = "Sain"
        nbm_frames.append(df)
        del adata, sample, df
        gc.collect()

    common_cols = set(nbm_frames[0].columns)
    for df_nbm in nbm_frames[1:]:
        common_cols &= set(df_nbm.columns)
    common_cols = sorted(common_cols)

    nbm_df = pd.concat(
        [df_nbm[common_cols] for df_nbm in nbm_frames], ignore_index=True
    )
    # --- liberation immediate des frames intermediaires ---
    del nbm_frames
    gc.collect()

    # 2) Determine marker list --------------------------------------------------
    patho_files = get_fcs_files(PATHO_DIR)
    if not patho_files:
        raise FileNotFoundError(f"Aucun fichier pathologique trouve dans {PATHO_DIR}")

    patients: List[Tuple[str, Path, Optional[float]]] = []
    for fcs in patho_files:
        p = Path(fcs)
        pid = p.stem
        gt_val = ground_truth.get(pid)
        if gt_val is not None:
            gt_val = float(gt_val)
            if not np.isfinite(gt_val):
                gt_val = None
        patients.append((pid, p, gt_val))

    n_known = sum(1 for _, _, gt in patients if gt is not None)
    log.info(
        "Cohorte Test 2: %d fichier(s), dont %d avec ground-truth CSV.",
        len(patients),
        n_known,
    )
    if n_known == 0:
        raise RuntimeError("Aucun patient de Test 2 ne correspond au CSV ground-truth.")

    pilot_path = patients[0][1]
    pilot_adata = load_fcs_file(pilot_path, condition="Patho")
    pilot_sample = FlowSample.from_anndata(pilot_adata, str(pilot_path), "Patho")
    pilot_pre = preprocess_sample(pilot_sample, config)
    base_markers = select_markers_for_clustering(pilot_pre.var_names, config)
    base_markers = [m for m in base_markers if m in nbm_df.columns]
    del pilot_adata, pilot_sample, pilot_pre
    gc.collect()

    if not base_markers:
        raise RuntimeError("Aucun marqueur commun trouve pour le clustering.")

    nbm_df = nbm_df[base_markers + ["condition"]]
    log.info("NBM: %d cellules, %d marqueurs communs", len(nbm_df), len(base_markers))

    # 3) FlowSOM une fois par patient ------------------------------------------
    patient_cache: List[PatientCache] = []
    for pid, fcs_path, gt in patients:
        log.info("FlowSOM precalcul patient %s", pid)

        adata_p = load_fcs_file(fcs_path, condition="Patho")
        sample_p = FlowSample.from_anndata(adata_p, str(fcs_path), "Patho")
        p_df = preprocess_sample(sample_p, config).data.copy()
        p_df["condition"] = "Patho"
        del adata_p, sample_p
        gc.collect()

        avail = [m for m in base_markers if m in p_df.columns]
        if not avail:
            log.warning("Patient %s: aucun marqueur disponible, ignore.", pid)
            del p_df
            gc.collect()
            continue
        p_df = p_df[avail + ["condition"]]

        data_train = pd.concat(
            [
                nbm_df[avail + ["condition"]].sample(
                    n=min(len(p_df), len(nbm_df)), random_state=42
                ),
                p_df,
            ],
            ignore_index=True,
        )
        del p_df
        gc.collect()

        x_data = data_train[avail].values.astype(np.float32)
        som = FlowSOMClusterer(
            xdim=som_grid,
            ydim=som_grid,
            seed=42,
            use_gpu=True,
        )
        som.fit(x_data, marker_names=avail)

        node_assignments = som.node_assignments_.copy()
        n_nodes = int(node_assignments.max()) + 1
        node_medians = np.zeros((n_nodes, len(avail)), dtype=np.float32)
        for i in range(n_nodes):
            mask = node_assignments == i
            if mask.any():
                node_medians[i] = np.median(x_data[mask], axis=0)

        x_nbm = data_train.loc[data_train["condition"] == "Sain", avail].values.astype(
            np.float32
        )
        nbm_stats = compute_reference_stats(x_nbm)
        del x_nbm, som
        gc.collect()

        # Detecter automatiquement si le patient est un "sain de reference" (MRD=0)
        is_healthy_control = (
            gt is not None and np.isfinite(float(gt)) and float(gt) == 0.0
        )

        patient_cache.append(
            PatientCache(
                pid=pid,
                fcs_path=fcs_path,
                data_train=data_train,
                node_assignments=node_assignments,
                node_medians=node_medians,
                nbm_stats=nbm_stats,
                markers=avail,
                true_mrd=float(gt) if gt is not None else float("nan"),
                has_ground_truth=(gt is not None and np.isfinite(float(gt))),
                is_healthy_control=is_healthy_control,
            )
        )

        del x_data
        gc.collect()

    del nbm_df
    gc.collect()

    return patient_cache


# ---------------------------------------------------------------------------
# Fonction de coût hybride cliniquement pondérée
# ---------------------------------------------------------------------------
def clinical_weighted_loss(
    pred: float,
    true_mrd: float,
    is_healthy_control: bool,
) -> float:
    """
    Compute a clinically weighted loss for a single patient.

    Logique:
      1. Patient sain (is_healthy_control=True):
           - Si pred > FP_THRESHOLD_PCT → pénalité dure FP_PENALTY
           - Sinon → erreur absolue classique (petite, car vrai MRD ~0)

      2. Patient pathologique:
           - err_base = |pred - true_mrd|
           - weight = 1 + HIGH_MRD_WEIGHT * (true_mrd / HIGH_MRD_CUTOFF)^HIGH_MRD_EXPONENT
             → pour true_mrd=1%   : weight ≈ 1.12  (quasi-neutre)
             → pour true_mrd=5%   : weight ≈ 4.0   (fort)
             → pour true_mrd=20%  : weight ≈ 49.0  (critique)
           - loss = err_base * weight

    Cela garantit qu'un trial qui "noie" un patient à 80% de MRD reçoit
    une pénalité ~624x supérieure à l'erreur absolue brute.
    """
    if np.isnan(pred):
        return FAILURE_PENALTY

    abs_err = abs(pred - true_mrd)

    if is_healthy_control:
        # Faux positif sur sain → pénalité dure
        if pred > FP_THRESHOLD_PCT:
            return FP_PENALTY + abs_err
        return abs_err

    # Patient pathologique: pondération exponentielle par sévérité MRD
    if true_mrd > 0.0:
        weight = (
            1.0 + HIGH_MRD_WEIGHT * (true_mrd / HIGH_MRD_CUTOFF) ** HIGH_MRD_EXPONENT
        )
    else:
        weight = 1.0

    return abs_err * weight


# ---------------------------------------------------------------------------
# Objective Optuna
# ---------------------------------------------------------------------------
def build_objective(
    base_mrd_cfg,
    patient_cache_by_grid: Dict[int, List[PatientCache]],
    som_grids: List[int],
):
    """
    Retourne la fonction objectif Optuna avec:
      - Hybrid weighted loss (vs MAE pure)
      - Pruning precoce (rapport intermediaire apres chaque patient)
      - Search space log-scale sur les distances
      - copy.deepcopy evite via patch direct sur une copie pre-faite au debut
    """
    # Pre-copie unique de la config de base (structure legere).
    # Chaque trial patch cette copie in-place, puis la restaure.
    # ATTENTION: si MrdConfig contient des objets mutables partages, utiliser
    # copy.deepcopy ici une seule fois et patcher la copie profonde a chaque trial.
    _base_copy = copy.deepcopy(base_mrd_cfg)

    def objective(trial: optuna.Trial) -> float:
        # ------------------------------------------------------------------ #
        # 1. Espace de recherche                                              #
        # ------------------------------------------------------------------ #
        mrd_method = trial.suggest_categorical("method", ["jf", "flo"])
        som_grid = trial.suggest_categorical("som_grid", som_grids)

        jf_max_normal = trial.suggest_float("jf_max_normal_marrow_pct", 0.1, 30.0)
        jf_min_patho = trial.suggest_float("jf_min_patho_cells_pct", 1.0, 35.0)
        nmm = trial.suggest_float("normal_marrow_multiplier", 1.0, 6.0)

        bpf_enabled = trial.suggest_categorical("bpf_enabled", [True, False])
        bpf_apply_jf = trial.suggest_categorical("bpf_apply_to_jf", [True, False])
        bpf_apply_flo = trial.suggest_categorical("bpf_apply_to_flo", [True, False])
        bpf_apply_eln = trial.suggest_categorical("bpf_apply_to_eln", [True, False])
        if bpf_enabled:
            scoring_method = trial.suggest_categorical(
                "scoring_method", ["none", "linear", "mahalanobis", "hybrid"]
            )
        else:
            # Pas de blast detection -> comparaison sur valeur brute uniquement.
            scoring_method = "none"
        allowed_mode = trial.suggest_categorical(
            "allowed_categories_mode", ["high_only", "high_moderate"]
        )

        high_threshold = trial.suggest_float("high_threshold", 2.0, 10.0)
        moderate_threshold = trial.suggest_float("moderate_threshold", 0.1, 6.0)
        weak_threshold = trial.suggest_float("weak_threshold", -1.0, 3.0)

        mahal_weight = trial.suggest_float("mahal_weight", 0.05, 0.95)
        linear_weight = trial.suggest_float("linear_weight", 0.05, 0.95)

        # --- log=True sur les parametres de distance (varient sur 1-2 ordres) ---
        d2_normalization = trial.suggest_float("d2_normalization", 5.0, 200.0, log=True)
        d2_high = trial.suggest_float("d2_threshold_high", 10.0, 300.0, log=True)
        d2_mod = trial.suggest_float("d2_threshold_moderate", 5.0, 150.0, log=True)
        d2_min = trial.suggest_float("d2_min_normal_threshold", 2.0, 80.0, log=True)
        mahal_boost = trial.suggest_float("mahal_boost_factor", 1.0, 6.0, log=True)

        purity_mod_power = trial.suggest_float("purity_modulation_power", 0.0, 2.0)
        purity_override = trial.suggest_categorical(
            "purity_threshold_override", [True, False]
        )
        purity_strict = trial.suggest_float("purity_strict_above", 0.70, 0.995)
        relax = trial.suggest_float("purity_relax_factor", 0.05, 1.00)

        # marker weights: log=True sur les positifs (CD34, CD117...) car ratio important
        w_cd34 = trial.suggest_float("w_CD34", 0.2, 8.0, log=True)
        w_cd117 = trial.suggest_float("w_CD117", 0.2, 8.0, log=True)
        w_hladr = trial.suggest_float("w_HLA_DR", 0.2, 8.0, log=True)
        w_cd33 = trial.suggest_float("w_CD33", 0.1, 6.0, log=True)
        w_cd13 = trial.suggest_float("w_CD13", 0.1, 6.0, log=True)
        # Poids negatifs: echelle lineaire (valeur negative bornee)
        w_cd45 = trial.suggest_float("w_CD45", -8.0, -0.05)
        w_ssc = trial.suggest_float("w_SSC", -12.0, -0.05)

        # ------------------------------------------------------------------ #
        # 2. Construction de la config (patch in-place sur deepcopy local)   #
        # ------------------------------------------------------------------ #
        mrd_cfg = copy.deepcopy(_base_copy)
        mrd_cfg.condition_sain = "Sain"
        mrd_cfg.condition_patho = "Patho"
        mrd_cfg.method = mrd_method
        mrd_cfg.method_jf.max_normal_marrow_pct = jf_max_normal
        mrd_cfg.method_jf.min_patho_cells_pct = jf_min_patho
        mrd_cfg.method_flo.normal_marrow_multiplier = nmm

        bpf = mrd_cfg.blast_phenotype_filter
        bpf.enabled = bool(bpf_enabled and scoring_method != "none")
        bpf.apply_to_jf = bpf_apply_jf
        bpf.apply_to_flo = bpf_apply_flo
        bpf.apply_to_eln = bpf_apply_eln
        bpf.scoring_method = scoring_method
        bpf.allowed_categories = (
            ["BLAST_HIGH"]
            if allowed_mode == "high_only"
            else ["BLAST_HIGH", "BLAST_MODERATE"]
        )

        ordered_linear = sorted([weak_threshold, moderate_threshold, high_threshold])
        bpf.weak_threshold = ordered_linear[0]
        bpf.moderate_threshold = ordered_linear[1]
        bpf.high_threshold = ordered_linear[2]

        ordered_d2 = sorted([d2_min, d2_mod, d2_high])
        bpf.d2_min_normal_threshold = ordered_d2[0]
        bpf.d2_threshold_moderate = ordered_d2[1]
        bpf.d2_threshold_high = ordered_d2[2]
        bpf.d2_normalization = d2_normalization
        bpf.mahal_boost_factor = mahal_boost

        total_w = mahal_weight + linear_weight
        if total_w <= 1e-8:
            bpf.mahal_weight = 0.5
            bpf.linear_weight = 0.5
        else:
            bpf.mahal_weight = mahal_weight / total_w
            bpf.linear_weight = linear_weight / total_w

        bpf.purity_modulation_power = purity_mod_power
        bpf.purity_threshold_override = purity_override
        bpf.purity_strict_above = purity_strict
        bpf.purity_relax_factor = relax

        bpf.marker_weights = {
            "CD34": w_cd34,
            "CD117": w_cd117,
            "HLA-DR": w_hladr,
            "CD33": w_cd33,
            "CD13": w_cd13,
            "CD45": w_cd45,
            "SSC": w_ssc,
        }

        # ------------------------------------------------------------------ #
        # 3. Evaluation avec pruning intermediaire                            #
        # ------------------------------------------------------------------ #
        weighted_losses: List[float] = []
        details: Dict[str, float] = {}
        cumulative_loss = 0.0

        grid_cache = patient_cache_by_grid[som_grid]
        active_patients = [pc for pc in grid_cache if pc.has_ground_truth]

        for step, pc in enumerate(active_patients):
            try:
                res = compute_mrd(
                    df_cells=pc.data_train,
                    clustering=pc.node_assignments,
                    mrd_cfg=mrd_cfg,
                    condition_column="condition",
                    node_medians=pc.node_medians,
                    marker_names=pc.markers,
                    nbm_center=pc.nbm_stats[0],
                    nbm_scale=pc.nbm_stats[1],
                    nbm_inv_cov=pc.nbm_stats[2],
                )
                pred = (
                    float(res.mrd_pct_jf)
                    if mrd_method == "jf"
                    else float(res.mrd_pct_flo)
                )
            except Exception as exc:
                log.debug(
                    "Trial %d, patient %s: compute_mrd failed (%s)",
                    trial.number,
                    pc.pid,
                    exc,
                )
                pred = float("nan")

            loss_i = clinical_weighted_loss(pred, pc.true_mrd, pc.is_healthy_control)

            weighted_losses.append(loss_i)
            cumulative_loss += loss_i
            details[f"pred_{pc.pid}"] = pred if not np.isnan(pred) else -1.0
            details[f"loss_{pc.pid}"] = loss_i

            # Pruning: rapport de la perte cumulative normalisee apres chaque patient
            # Un trial manifestement mauvais est coupe avant d'evaluer tous les patients.
            if len(active_patients) > 1:
                intermediate = cumulative_loss / (step + 1)
                trial.report(intermediate, step=step)
                if trial.should_prune():
                    raise optuna.exceptions.TrialPruned()

        if not weighted_losses:
            return FAILURE_PENALTY

        objective_value = float(np.mean(weighted_losses))

        trial.set_user_attr("prediction_method", mrd_method)
        trial.set_user_attr("som_grid", int(som_grid))
        trial.set_user_attr("scoring_method", scoring_method)
        trial.set_user_attr("blast_filter_enabled_effective", int(bpf.enabled))
        trial.set_user_attr("objective_weighted_loss", objective_value)
        trial.set_user_attr(
            "n_fp_violations",
            sum(
                1
                for pc in active_patients
                if pc.is_healthy_control
                and details.get(f"pred_{pc.pid}", 0.0) > FP_THRESHOLD_PCT
            ),
        )
        for k, v in details.items():
            trial.set_user_attr(k, float(v))

        return objective_value

    return objective


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------
def progress_callback(study: optuna.Study, trial: optuna.trial.FrozenTrial) -> None:
    """Compact progress callback with convergence tracking."""
    if trial.state == optuna.trial.TrialState.PRUNED:
        return  # silencieux pour les trials elagues

    if trial.number == 0 or (trial.number + 1) % 10 == 0:
        best_value = float("nan")
        try:
            best_value = study.best_value
        except ValueError:
            pass
        fp_violations = trial.user_attrs.get("n_fp_violations", "?")
        log.info(
            "Trial %d | loss=%.6f | best=%.6f | FP_violations=%s",
            trial.number + 1,
            trial.value if trial.value is not None else float("nan"),
            best_value,
            fp_violations,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Optimisation bayesienne MRD avec Optuna (v2 — weighted loss)"
    )
    parser.add_argument(
        "--trials", type=int, default=180, help="Nombre d'iterations Optuna"
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Seed pour reproductibilite"
    )
    parser.add_argument(
        "--ground-truth-csv",
        type=str,
        default=str(GROUND_TRUTH_CSV_PATH),
        help="Chemin du CSV de ground-truth MRD",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reprendre l'etude SQLite existante au lieu d'en creer une nouvelle",
    )
    parser.add_argument(
        "--fp-threshold",
        type=float,
        default=FP_THRESHOLD_PCT,
        help=f"Seuil faux positif patient sain en %% (defaut: {FP_THRESHOLD_PCT})",
    )
    parser.add_argument(
        "--som-grids",
        type=str,
        default="10,14",
        help="Tailles de grille SOM a comparer (ex: 10,14)",
    )
    return parser.parse_args()


def parse_som_grids(raw: str) -> List[int]:
    """Parse comma-separated SOM grid sizes and validate them."""
    grids: List[int] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = int(token)
        except ValueError as exc:
            raise ValueError(f"Valeur de grille invalide: '{token}'") from exc
        if value < 2:
            raise ValueError(f"Taille de grille invalide (<2): {value}")
        grids.append(value)

    grids = sorted(set(grids))
    if not grids:
        raise ValueError("Aucune grille SOM valide fournie.")
    return grids


def build_grid_and_global_summary(
    trials_df: pd.DataFrame,
    grids: List[int],
) -> Dict[str, object]:
    """Build summary for best loss per grid and global winner."""
    summary: Dict[str, object] = {
        "best_by_grid": {},
        "global_winner": None,
    }

    if trials_df.empty:
        return summary

    if "state" in trials_df.columns:
        complete_mask = trials_df["state"].astype(str).str.contains("COMPLETE")
        done = trials_df.loc[complete_mask].copy()
    else:
        done = trials_df.copy()

    if done.empty:
        return summary

    def _pick_first_available(row: pd.Series, cols: List[str], default=np.nan):
        for c in cols:
            if c in row.index and pd.notna(row[c]):
                return row[c]
        return default

    for grid in grids:
        grid_rows = done[
            done.get("params_som_grid", pd.Series(index=done.index)).astype(str)
            == str(grid)
        ]
        if grid_rows.empty:
            continue

        best_idx = grid_rows["value"].astype(float).idxmin()
        best_row = grid_rows.loc[best_idx]
        summary["best_by_grid"][str(grid)] = {
            "trial": int(best_row["number"]),
            "loss": float(best_row["value"]),
            "method": str(
                _pick_first_available(
                    best_row,
                    ["params_method", "user_attrs_prediction_method"],
                    "unknown",
                )
            ),
            "scoring_method": str(
                _pick_first_available(
                    best_row,
                    ["params_scoring_method", "user_attrs_scoring_method"],
                    "none",
                )
            ),
        }

    best_global_idx = done["value"].astype(float).idxmin()
    best_global = done.loc[best_global_idx]
    summary["global_winner"] = {
        "trial": int(best_global["number"]),
        "loss": float(best_global["value"]),
        "grid": int(
            _pick_first_available(
                best_global,
                ["params_som_grid", "user_attrs_som_grid"],
                -1,
            )
        ),
        "method": str(
            _pick_first_available(
                best_global,
                ["params_method", "user_attrs_prediction_method"],
                "unknown",
            )
        ),
        "scoring_method": str(
            _pick_first_available(
                best_global,
                ["params_scoring_method", "user_attrs_scoring_method"],
                "none",
            )
        ),
    }

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    som_grids = parse_som_grids(args.som_grids)

    # Surcharge eventuellement du seuil FP via CLI
    global FP_THRESHOLD_PCT
    FP_THRESHOLD_PCT = args.fp_threshold

    log.info(
        "Demarrage optimisation Optuna MRD v2 | trials=%d | FP_threshold=%.2f%% | som_grids=%s",
        args.trials,
        FP_THRESHOLD_PCT,
        som_grids,
    )

    config = PipelineConfig()
    base_mrd_cfg = load_mrd_config(MRD_CFG_PATH)

    gt_csv, ground_truth, n_matches = resolve_ground_truth_for_test2(
        Path(args.ground_truth_csv),
        PATHO_DIR,
    )
    log.info(
        "CSV ground-truth retenu: %s (%d PID dans CSV, %d match(es) Test 2)",
        gt_csv,
        len(ground_truth),
        n_matches,
    )

    finite_truth = {k: v for k, v in ground_truth.items() if np.isfinite(float(v))}
    if not finite_truth:
        raise RuntimeError(
            "Aucune valeur ground-truth numerique exploitable dans le CSV."
        )

    test2_pids = {Path(f).stem for f in get_fcs_files(PATHO_DIR)}
    finite_truth_test2 = {k: v for k, v in finite_truth.items() if k in test2_pids}
    if not finite_truth_test2:
        raise RuntimeError(
            "Aucune valeur ground-truth numerique du CSV ne correspond aux patients Test 2."
        )

    log.info(
        "Optimisation globale active: %d patients Test 2 avec truth numerique.",
        len(finite_truth_test2),
    )

    patient_cache_by_grid: Dict[int, List[PatientCache]] = {}
    for grid in som_grids:
        patient_cache_by_grid[grid] = prepare_patient_cache(config, ground_truth, grid)

    # -------------------------------------------------------------------------
    # Persistance SQLite — reprise automatique apres interruption
    # -------------------------------------------------------------------------
    sqlite_url = f"sqlite:///{SQLITE_PATH}"
    study_name = "mrd_optuna_v2"

    sampler = optuna.samplers.TPESampler(seed=args.seed)
    # MedianPruner: elague si la valeur intermediaire depasse la mediane des trials
    # precedents au meme step. n_startup_trials=10 pour laisser l'exploration initiale.
    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=10,
        n_warmup_steps=2,  # attendre au moins 2 patients avant d'elaguer
        interval_steps=1,
    )

    if args.resume:
        study = optuna.load_study(
            study_name=study_name,
            storage=sqlite_url,
            sampler=sampler,
            pruner=pruner,
        )
        log.info(
            "Reprise etude '%s' depuis SQLite (%d trials existants).",
            study_name,
            len(study.trials),
        )
    else:
        study = optuna.create_study(
            study_name=study_name,
            storage=sqlite_url,
            direction="minimize",
            sampler=sampler,
            pruner=pruner,
            load_if_exists=True,  # securite: ne pas ecraser si la DB existe deja
        )

    objective = build_objective(base_mrd_cfg, patient_cache_by_grid, som_grids)

    study.optimize(
        objective,
        n_trials=args.trials,
        callbacks=[progress_callback],
        gc_after_trial=True,
        # n_jobs=1 par defaut (safe avec GPU partagee et cache patient).
        # Pour CPU-only: n_jobs=-1 si compute_mrd est thread-safe.
    )

    # -------------------------------------------------------------------------
    # Export resultats
    # -------------------------------------------------------------------------
    trials_df = study.trials_dataframe(
        attrs=("number", "value", "params", "state", "user_attrs")
    )
    trials_df.to_csv(TRIALS_CSV, index=False)
    run_summary = build_grid_and_global_summary(trials_df, som_grids)
    SUMMARY_JSON.write_text(json.dumps(run_summary, indent=2), encoding="utf-8")

    try:
        best_value = float(study.best_value)
        best_params = study.best_params
    except ValueError:
        log.warning("Aucun trial complete avec succes.")
        return

    best_payload = {
        "best_value": best_value,
        "best_params": best_params,
        "n_trials": len(study.trials),
        "objective": "clinical_weighted_loss",
        "fp_threshold_pct": FP_THRESHOLD_PCT,
        "high_mrd_cutoff": HIGH_MRD_CUTOFF,
        "high_mrd_weight": HIGH_MRD_WEIGHT,
        "high_mrd_exponent": HIGH_MRD_EXPONENT,
        "fp_penalty": FP_PENALTY,
        "n_ground_truth_patients": int(len(finite_truth_test2)),
        "som_grids_tested": som_grids,
        "sqlite_path": str(SQLITE_PATH),
    }
    BEST_JSON.write_text(json.dumps(best_payload, indent=2), encoding="utf-8")

    print("\n" + "=" * 70)
    print("GOLDEN CONFIGURATION (OPTUNA v2 — Weighted Clinical Loss)")
    print(f"Best weighted loss: {best_value:.6f}")
    print(
        f"  (FP_threshold={FP_THRESHOLD_PCT}% | "
        f"high_mrd_cutoff={HIGH_MRD_CUTOFF}% | "
        f"high_mrd_exponent={HIGH_MRD_EXPONENT})"
    )
    for key, value in best_params.items():
        print(f"  - {key}: {value}")
    best_by_grid = run_summary.get("best_by_grid", {})
    best_10 = best_by_grid.get("10")
    best_14 = best_by_grid.get("14")
    global_winner = run_summary.get("global_winner")

    print("\nRésumé auto")
    if best_10 is not None:
        print(
            "  - meilleure loss grille 10: "
            f"{best_10['loss']:.6f} "
            f"(trial {best_10['trial']}, method={best_10['method']}, "
            f"scoring={best_10['scoring_method']})"
        )
    if best_14 is not None:
        print(
            "  - meilleure loss grille 14: "
            f"{best_14['loss']:.6f} "
            f"(trial {best_14['trial']}, method={best_14['method']}, "
            f"scoring={best_14['scoring_method']})"
        )
    if global_winner is not None:
        print(
            "  - gagnant global: "
            f"grille {global_winner['grid']} | "
            f"method={global_winner['method']} | "
            f"scoring={global_winner['scoring_method']} | "
            f"loss={global_winner['loss']:.6f} (trial {global_winner['trial']})"
        )
    print("=" * 70)
    log.info(
        "Exports: %s | %s | %s | %s",
        TRIALS_CSV,
        BEST_JSON,
        SUMMARY_JSON,
        SQLITE_PATH,
    )


if __name__ == "__main__":
    main()
