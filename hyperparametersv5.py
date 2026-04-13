"""
hyperparametersv5.py — Optimisation Clinique Finale
===================================================
Correction du ModuleNotFoundError et optimisation du scoring hybride.
"""

import os
import sys
import copy
import gc
import logging
import time
import numpy as np
import pandas as pd
from itertools import product
from pathlib import Path
from dataclasses import dataclass

# ── FIX DU PYTHONPATH (À placer avant les imports src) ──────────────────────
# CURRENT_DIR = flowsom_pipeline_pro/
# PARENT_DIR  = Perplexity/  ← nécessaire pour résoudre "flowsom_pipeline_pro.xxx"
CURRENT_DIR = Path(__file__).resolve().parent
PARENT_DIR = CURRENT_DIR.parent

if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

os.environ["PYTHONPATH"] = str(PARENT_DIR)

# ── Imports internes (Maintenant sécurisés) ──────────────────────────────────
try:
    from src.io.fcs_reader import load_fcs_file, get_fcs_files
    from src.models.sample import FlowSample
    from src.services.preprocessing_service import preprocess_sample
    from src.services.clustering_service import select_markers_for_clustering
    from src.core.clustering import FlowSOMClusterer
    from src.analysis.mrd_calculator import load_mrd_config, compute_mrd
    from src.analysis.blast_detection import compute_reference_stats
    from config.pipeline_config import PipelineConfig
except ImportError as e:
    print(f"\n❌ Erreur d'importation : {e}")
    print(f"DEBUG - sys.path : {sys.path}")
    sys.exit(1)

# ── Configuration des Chemins ────────────────────────────────────────────────
HEALTHY_DIR = Path("F:/BLAST110/BLAST110/FCS/NBM/T1")
PATHO_DIR = Path("F:/BLAST110/BLAST110/FCS/Test")
OUTPUT_DIR = CURRENT_DIR / "Résults" / "Optimisation_V5"
OUTPUT_CSV = OUTPUT_DIR / "results_grid_search_v5.csv"
MRD_CFG_PATH = CURRENT_DIR / "config" / "mrd_config.yaml"

# ── Ground Truth (Vérité Terrain) ───────────────────────────────────────────
GROUND_TRUTH = {
    "BLAST110_105_P1": 0.27,  # Sain
    "BLAST110_25_P1": 82.56,  # Malade atypique
    "BLAST110_57_P1": 91.02,  # Malade atypique
}

# ── Grille de Recherche (Grid Search) ───────────────────────────────────────
PARAM_GRID = {
    "d2_threshold_moderate": [20.0, 25.0, 30.0],
    "d2_min_normal_threshold": [10.0, 12.0, 15.0],
    "mahal_boost_factor": [1.5, 2.0],
    "purity_relax_factor": [0.7, 0.8, 0.9],
    "normal_marrow_multiplier": [2.0, 3.0],
}

# ── Paramètres Fixes (Réglage Chirurgical) ──────────────────────────────────
MARKER_WEIGHTS = {
    "CD34": 3.0,
    "CD117": 2.5,
    "HLA-DR": 2.0,
    "CD33": 1.5,
    "CD13": 1.0,
    "CD45": -2.0,
    "SSC": -5.0,
}

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


@dataclass
class PatientCache:
    pid: str
    data_train: pd.DataFrame
    node_assignments: np.ndarray
    node_medians: np.ndarray
    nbm_stats: tuple
    true_mrd: float


def run_v5_optimization():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config = PipelineConfig()
    base_mrd_cfg = load_mrd_config(MRD_CFG_PATH)

    # 1. Charger la Moelle Saine (NBM) — colonnes communes (pattern v4)
    nbm_files = get_fcs_files(HEALTHY_DIR)
    nbm_frames = []
    log.info("=== Préparation du pool NBM de référence ===")
    for f in nbm_files[:10]:
        adata = load_fcs_file(f, condition="Sain")
        s = FlowSample.from_anndata(adata, str(f), "Sain")
        df = preprocess_sample(s, config).data.copy()
        df["condition"] = "Sain"
        nbm_frames.append(df)
        del s, adata
        gc.collect()

    # Colonnes communes entre tous les fichiers NBM (évite KeyError sur colonnes absentes)
    common_cols = set(nbm_frames[0].columns)
    for df in nbm_frames[1:]:
        common_cols &= set(df.columns)
    common_cols = sorted(common_cols)
    nbm_df = pd.concat([df[common_cols] for df in nbm_frames], ignore_index=True)
    del nbm_frames
    gc.collect()
    log.info("  NBM total : %d cellules, colonnes : %s", len(nbm_df), list(common_cols))

    # Marqueurs déterminés dynamiquement depuis le premier patient (pattern v4)
    pilot_path = PATHO_DIR / f"{next(iter(GROUND_TRUTH))}.fcs"
    pilot_adata = load_fcs_file(pilot_path, condition="Patho")
    pilot_sample = FlowSample.from_anndata(pilot_adata, str(pilot_path), "Patho")
    pilot_preprocessed = preprocess_sample(pilot_sample, config)
    markers = select_markers_for_clustering(pilot_preprocessed.var_names, config)
    markers = [m for m in markers if m in nbm_df.columns]
    log.info("Marqueurs (%d) : %s", len(markers), markers)
    del pilot_adata, pilot_sample, pilot_preprocessed
    gc.collect()

    # Restreindre nbm_df aux marqueurs utiles + condition
    nbm_df = nbm_df[[m for m in markers if m in nbm_df.columns] + ["condition"]]

    # 2. PASSE 1 : Précalcul FlowSOM pour les 3 patients critiques
    patient_cache = []
    for pid, true_val in GROUND_TRUTH.items():
        log.info(f"Passe 1 : Clustering FlowSOM pour {pid}...")
        fcs_path = PATHO_DIR / f"{pid}.fcs"
        adata_p = load_fcs_file(fcs_path, condition="Patho")
        p_sample = FlowSample.from_anndata(adata_p, str(fcs_path), "Patho")
        p_df = preprocess_sample(p_sample, config).data.copy()
        p_df["condition"] = "Patho"
        del adata_p, p_sample
        gc.collect()

        # Filtrer aux marqueurs disponibles dans ce patient
        avail = [m for m in markers if m in p_df.columns]
        p_df = p_df[avail + ["condition"]]

        # Ratio 1:1 pour Amsterdam style
        data_train = pd.concat(
            [nbm_df.sample(n=min(len(p_df), len(nbm_df)), random_state=42), p_df],
            ignore_index=True,
        )

        X = data_train[avail].values.astype(np.float32)
        som = FlowSOMClusterer(xdim=14, ydim=14, seed=42, use_gpu=True)
        som.fit(X, marker_names=avail)
        node_assignments = som.node_assignments_.copy()

        n_nodes = int(node_assignments.max()) + 1
        node_medians = np.zeros((n_nodes, len(avail)), dtype=np.float32)
        for i in range(n_nodes):
            mask = node_assignments == i
            if mask.any():
                node_medians[i] = np.median(X[mask], axis=0)

        X_nbm = data_train.loc[data_train["condition"] == "Sain", avail].values.astype(np.float32)
        nbm_stats = compute_reference_stats(X_nbm)
        del X, som
        gc.collect()

        patient_cache.append(
            PatientCache(
                pid=pid,
                data_train=data_train,
                node_assignments=node_assignments,
                node_medians=node_medians,
                nbm_stats=nbm_stats,
                true_mrd=true_val,
            )
        )

    # 3. PASSE 2 : Grid Search
    combos = list(product(*PARAM_GRID.values()))
    results = []
    log.info(f"Passe 2 : Test de {len(combos)} combinaisons d'hyperparamètres...")

    for d2_mod, d2_min, boost, relax, nmm in combos:
        mrd_cfg = copy.deepcopy(base_mrd_cfg)
        mrd_cfg.blast_phenotype_filter.scoring_method = "hybrid"
        mrd_cfg.blast_phenotype_filter.d2_threshold_moderate = d2_mod
        mrd_cfg.blast_phenotype_filter.d2_min_normal_threshold = d2_min
        mrd_cfg.blast_phenotype_filter.mahal_boost_factor = boost
        mrd_cfg.blast_phenotype_filter.purity_relax_factor = relax
        mrd_cfg.method_flo.normal_marrow_multiplier = nmm
        mrd_cfg.blast_phenotype_filter.marker_weights = MARKER_WEIGHTS

        for pc in patient_cache:
            avail = [c for c in pc.data_train.columns if c != "condition"]
            try:
                res = compute_mrd(
                    df_cells=pc.data_train,
                    clustering=pc.node_assignments,
                    mrd_cfg=mrd_cfg,
                    condition_column="condition",
                    node_medians=pc.node_medians,
                    marker_names=avail,
                    nbm_center=pc.nbm_stats[0],
                    nbm_scale=pc.nbm_stats[1],
                    nbm_inv_cov=pc.nbm_stats[2],
                )
                mrd_pred = res.mrd_pct_flo
            except Exception as exc:
                log.error("  ✗ pid=%s : %s", pc.pid, exc)
                mrd_pred = float("nan")

            error = abs(mrd_pred - pc.true_mrd) if not (mrd_pred != mrd_pred) else float("nan")
            clin_score = error * 3.0 if pc.true_mrd < 1.0 else error

            results.append(
                {
                    "pid": pc.pid,
                    "d2_mod": d2_mod,
                    "d2_min": d2_min,
                    "boost": boost,
                    "relax": relax,
                    "nmm": nmm,
                    "predicted": mrd_pred,
                    "true": pc.true_mrd,
                    "error": error,
                    "clinical_score": clin_score,
                }
            )

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)

    summary = (
        df.groupby(["d2_mod", "d2_min", "boost", "relax", "nmm"])["clinical_score"]
        .mean()
        .reset_index()
        .sort_values("clinical_score")
    )
    print("\n=== TOP 5 DES CONFIGURATIONS CLINIQUES (V5) ===")
    print(summary.head(5).to_string(index=False))


if __name__ == "__main__":
    run_v5_optimization()
