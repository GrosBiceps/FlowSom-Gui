"""
flowsom_plots.py — Visualisations FlowSOM (MST, heatmap MFI, UMAP).

Génère les graphiques spécifiques à l'analyse FlowSOM:
  - Heatmap des profils d'expression par métacluster
  - UMAP coloré par métacluster (attention au réglage de max_pts pour les gros datasets)
  - Graphique en barres de la taille des métaclusters
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.figure as _mpl_figure

    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False

try:
    import seaborn as sns

    _SNS_AVAILABLE = True
except ImportError:
    _SNS_AVAILABLE = False

from flowsom_pipeline_pro.src.visualization.plot_helpers import (
    BG_COLOR,
    TEXT_COLOR,
    SPINE_COLOR,
    apply_dark_style,
    save_figure,
)
from flowsom_pipeline_pro.src.utils.logger import get_logger

_logger = get_logger("visualization.flowsom_plots")


def plot_mfi_heatmap(
    mfi_matrix: pd.DataFrame,
    output_path: Path | str,
    *,
    figsize: Tuple[int, int] = (16, 8),
    cmap: str = "magma",
    title: str = "Profils d'expression par métacluster (MFI normalisée)",
) -> Optional["matplotlib.figure.Figure"]:
    """
    Heatmap des profils d'expression MFI par métacluster.

    Chaque ligne = un métacluster, chaque colonne = un marqueur.
    Les valeurs sont normalisées par marqueur (z-score) pour
    rendre toutes les intensités comparables visuellement.

    Args:
        mfi_matrix: DataFrame [n_clusters × n_markers] avec les MFI.
        output_path: Chemin PNG de sortie.
        figsize: Taille de la figure.
        cmap: Colormap matplotlib.
        title: Titre du graphique.

    Returns:
        True si succès.
    """
    if not _MPL_AVAILABLE or not _SNS_AVAILABLE:
        _logger.warning("matplotlib et seaborn requis pour plot_mfi_heatmap")
        return None

    try:
        # Normalisation par marqueur (z-score colonne)
        mfi_norm = mfi_matrix.copy().astype(float)
        stds = mfi_norm.std()
        stds[stds == 0] = 1.0
        mfi_norm = (mfi_norm - mfi_norm.mean()) / stds

        fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)

        sns.heatmap(
            mfi_norm,
            ax=ax,
            cmap=cmap,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            linecolor=SPINE_COLOR,
            cbar_kws={"label": "z-score MFI", "shrink": 0.8},
        )

        ax.set_title(title, fontsize=14, fontweight="bold", color=TEXT_COLOR, pad=15)
        ax.set_xlabel("Marqueurs", fontsize=12, color=TEXT_COLOR)
        ax.set_ylabel("Métacluster", fontsize=12, color=TEXT_COLOR)
        ax.tick_params(colors=TEXT_COLOR, labelsize=10, rotation=45)
        ax.set_facecolor(BG_COLOR)

        # Style colorbar
        cbar = ax.collections[0].colorbar
        if cbar:
            cbar.ax.tick_params(colors=TEXT_COLOR)
            cbar.ax.yaxis.label.set_color(TEXT_COLOR)

        save_figure(fig, output_path)
        _logger.info("Heatmap MFI sauvegardée: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_mfi_heatmap: %s", exc)
        return None


def plot_metacluster_sizes(
    metaclustering: np.ndarray,
    n_metaclusters: int,
    output_path: Path | str,
    condition_labels: Optional[np.ndarray] = None,
    title: str = "Distribution des cellules par métacluster",
) -> Optional["matplotlib.figure.Figure"]:
    """
    Figure 1×2 : pie chart (gauche) + bar chart horizontal (droite).

    Fidèle à la section §15 du pipeline monolithique flowsom_pipeline.py :
    le panneau gauche est un camembert des proportions, le panneau droit
    un bar chart horizontal des tailles absolues.
    Si condition_labels est fourni, le panneau droit utilise des barres
    empilées par condition (Sain/Patho).

    Args:
        metaclustering: Assignation par cellule (n_cells,).
        n_metaclusters: Nombre total de métaclusters.
        output_path: Chemin PNG de sortie.
        condition_labels: Labels de condition par cellule (optionnel).
        title: Titre du graphique.

    Returns:
        Figure matplotlib ou None si matplotlib absent.
    """
    if not _MPL_AVAILABLE:
        return None

    try:
        cluster_ids = np.arange(n_metaclusters)
        counts_total = np.array([int((metaclustering == i).sum()) for i in cluster_ids])
        labels = [f"MC{i}" for i in cluster_ids]
        colors = plt.cm.tab20(np.linspace(0, 1, n_metaclusters))

        fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor=BG_COLOR)

        # ── Panneau gauche : Pie chart ───────────────────────────────────────
        ax_pie = axes[0]
        ax_pie.set_facecolor(BG_COLOR)
        wedges, texts, autotexts = ax_pie.pie(
            counts_total,
            labels=labels,
            colors=colors,
            autopct="%1.1f%%",
            pctdistance=0.8,
            textprops={"color": TEXT_COLOR, "fontsize": 9},
        )
        for at in autotexts:
            at.set_color(TEXT_COLOR)
            at.set_fontsize(8)
        ax_pie.set_title(
            "Distribution des Cellules par Métacluster",
            fontsize=12,
            fontweight="bold",
            color=TEXT_COLOR,
        )

        # ── Panneau droite : Bar chart ───────────────────────────────────────
        ax_bar = axes[1]
        ax_bar.set_facecolor(BG_COLOR)

        if condition_labels is None:
            ax_bar.barh(
                range(n_metaclusters), counts_total, color=colors, edgecolor="white"
            )
            ax_bar.set_xlabel("Nombre de cellules", fontsize=11, color=TEXT_COLOR)
        else:
            unique_conds = sorted(set(condition_labels))
            cond_colors = ["#89b4fa", "#f38ba8", "#a6e3a1", "#f9e2af"]
            left = np.zeros(n_metaclusters)
            for ci, cond in enumerate(unique_conds):
                cond_mask = condition_labels == cond
                cnts = np.array(
                    [
                        int(((metaclustering == k) & cond_mask).sum())
                        for k in cluster_ids
                    ]
                )
                ax_bar.barh(
                    range(n_metaclusters),
                    cnts,
                    left=left,
                    color=cond_colors[ci % len(cond_colors)],
                    edgecolor="white",
                    linewidth=0.5,
                    label=str(cond),
                )
                left += cnts
            ax_bar.legend(
                facecolor="#313244",
                labelcolor=TEXT_COLOR,
                edgecolor=SPINE_COLOR,
                fontsize=9,
            )
            ax_bar.set_xlabel("Nombre de cellules", fontsize=11, color=TEXT_COLOR)

        ax_bar.set_yticks(range(n_metaclusters))
        ax_bar.set_yticklabels(labels, fontsize=9, color=TEXT_COLOR)
        ax_bar.set_title(
            "Taille des Métaclusters", fontsize=12, fontweight="bold", color=TEXT_COLOR
        )
        ax_bar.grid(axis="x", alpha=0.3, linestyle="--", color=SPINE_COLOR)
        ax_bar.tick_params(colors=TEXT_COLOR)
        for sp in ax_bar.spines.values():
            sp.set_color(SPINE_COLOR)

        plt.tight_layout()
        save_figure(fig, output_path)
        _logger.info("Distribution métaclusters sauvegardée: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_metacluster_sizes: %s", exc)
        return None


def plot_umap(
    umap_coords: np.ndarray,
    metaclustering: np.ndarray,
    output_path: Path | str,
    n_metaclusters: int = 10,
    title: str = "UMAP — Coloré par métacluster FlowSOM",
    max_pts: int = 10000,
    seed: int = 42,
) -> Optional["matplotlib.figure.Figure"]:
    """
    Scatter UMAP coloré par métacluster.

    Args:
        umap_coords: Coordonnées UMAP (n_cells, 2).
        metaclustering: Assignation de métacluster (n_cells,).
        output_path: Chemin PNG de sortie.
        n_metaclusters: Nombre de métaclusters (pour colormap discrète).
        title: Titre du graphique.
        max_pts: Sous-échantillonnage si > max_pts cellules.
        seed: Graine pour sous-échantillonnage.

    Returns:
        True si succès.
    """
    if not _MPL_AVAILABLE:
        return None

    try:
        if len(umap_coords) > max_pts:
            rng = np.random.default_rng(seed)
            idx = rng.choice(len(umap_coords), max_pts, replace=False)
            umap_coords = umap_coords[idx]
            metaclustering = metaclustering[idx]

        # Colormap discrète Tab20 (jusqu'à 20 couleurs distinctes)
        cmap = plt.cm.get_cmap("tab20", n_metaclusters)
        colors = [cmap(i % 20) for i in metaclustering]

        fig, ax = plt.subplots(figsize=(12, 10), facecolor=BG_COLOR)

        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=metaclustering,
            cmap="tab20",
            vmin=0,
            vmax=n_metaclusters,
            s=3,
            alpha=0.5,
            edgecolors="none",
            rasterized=True,
        )

        ax.set_xlabel("UMAP 1", fontsize=13, fontweight="bold", color=TEXT_COLOR)
        ax.set_ylabel("UMAP 2", fontsize=13, fontweight="bold", color=TEXT_COLOR)
        ax.set_title(title, fontsize=14, fontweight="bold", color=TEXT_COLOR, pad=12)
        ax.set_facecolor(BG_COLOR)
        ax.tick_params(colors=TEXT_COLOR)

        for spine in ax.spines.values():
            spine.set_color(SPINE_COLOR)

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label("Métacluster", color=TEXT_COLOR, fontsize=11)
        cbar.ax.tick_params(colors=TEXT_COLOR)
        cbar.set_ticks(np.arange(n_metaclusters))
        cbar.set_ticklabels([f"MC{i}" for i in range(n_metaclusters)])

        save_figure(fig, output_path)
        _logger.info("UMAP sauvegardé: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_umap: %s", exc)
        return None


# ─────────────────────────────────────────────────────────────────────────────
#  Utilitaires de jitter pour les graphiques FlowSOM-style
# ─────────────────────────────────────────────────────────────────────────────


def circular_jitter(
    n_points: int,
    cluster_ids: np.ndarray,
    node_sizes: np.ndarray,
    max_radius: float = 0.45,
    min_radius: float = 0.10,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Jitter circulaire vectorisé pour représentation FlowSOM-R.

    Chaque cellule est décalée d'une distance aléatoire dans un disque
    centré sur son nœud SOM. Le rayon maximum du disque est proportionnel
    à la racine carrée de la taille du nœud (conforme à FlowSOM R package).

    Formule du rayon:
        r_max = min_radius + (max_radius - min_radius) × √(node_size / max_node_size)

    Args:
        n_points: Nombre de cellules.
        cluster_ids: Assignation de nœud par cellule (n_cells,), entiers.
        node_sizes: Taille (n_cells) de chaque nœud SOM (n_nodes,).
        max_radius: Rayon maximum pour le nœud le plus grand (défaut 0.45).
        min_radius: Rayon minimum pour les petits nœuds (défaut 0.10).
        seed: Graine aléatoire pour la reproductibilité.

    Returns:
        Tuple (jitter_x, jitter_y) de dtype float32.
    """
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0.0, 2.0 * np.pi, n_points)
    u = rng.uniform(0.0, 1.0, n_points)

    ids = cluster_ids.astype(int)
    max_size_val = float(node_sizes.max())
    if max_size_val <= 0.0:
        max_size_val = 1.0

    radii = min_radius + (max_radius - min_radius) * np.sqrt(
        node_sizes[ids] / max_size_val
    )
    # √u pour une densité uniforme dans le disque (pas une concentration au centre)
    r = np.sqrt(u) * radii

    jitter_x = (r * np.cos(theta)).astype(np.float32)
    jitter_y = (r * np.sin(theta)).astype(np.float32)
    return jitter_x, jitter_y


# ─────────────────────────────────────────────────────────────────────────────
#  MST matplotlib statique
# ─────────────────────────────────────────────────────────────────────────────


def plot_mst_static(
    clusterer: Any,
    mfi_matrix: pd.DataFrame,
    metaclustering: np.ndarray,
    output_path: Path | str,
    *,
    figsize: Tuple[int, int] = (14, 12),
    title: str = "Arbre MST FlowSOM — Topologie des nodes SOM",
) -> Optional["matplotlib.figure.Figure"]:
    """
    Graphique MST statique matplotlib.

    Nodes colorés par métacluster dominant, taille proportionnelle au nombre
    de cellules. Arêtes de l'arbre minimal (MST sur le codebook).

    Args:
        clusterer: FlowSOMClusterer après fit().
        mfi_matrix: DataFrame [n_metaclusters × n_markers] — MFI médiane.
        metaclustering: Assignation métacluster par cellule (n_cells,).
        output_path: Chemin PNG de sortie.
        figsize: Taille de la figure.
        title: Titre du graphique.

    Returns:
        Figure matplotlib ou None si échec.
    """
    if not _MPL_AVAILABLE:
        return None

    try:
        from scipy.spatial.distance import cdist
        from scipy.sparse.csgraph import minimum_spanning_tree
        from scipy.sparse import csr_matrix
        from collections import defaultdict

        layout_coords = clusterer.get_layout_coords()  # (n_nodes, 2)
        node_sizes = clusterer.get_node_sizes()  # (n_nodes,)
        n_nodes = clusterer.n_nodes

        # ── Métacluster dominant par node ─────────────────────────────────────
        mc_per_node = np.full(n_nodes, -1, dtype=int)
        na = getattr(clusterer, "node_assignments_", None)
        ma = getattr(clusterer, "metacluster_assignments_", None)
        if na is not None and ma is not None:
            for node_id in range(n_nodes):
                cells_in_node = na == node_id
                if cells_in_node.any():
                    mc_per_node[node_id] = int(np.bincount(ma[cells_in_node]).argmax())

        # ── Arêtes MST depuis le codebook ─────────────────────────────────────
        edges: List[Tuple[int, int]] = []
        codebook = None
        fsm = getattr(clusterer, "_fsom_model", None)
        if fsm is not None:
            if hasattr(fsm, "codes"):
                codebook = np.asarray(fsm.codes, dtype=float)
            elif hasattr(fsm, "model") and hasattr(fsm.model, "codes"):
                codebook = np.asarray(fsm.model.codes, dtype=float)
        if codebook is not None:
            dist_mat = cdist(codebook, codebook, metric="euclidean")
            mst = minimum_spanning_tree(csr_matrix(dist_mat))
            coo = mst.tocoo()
            edges = list(zip(coo.row.tolist(), coo.col.tolist()))

        # ── Dessin ────────────────────────────────────────────────────────────
        n_meta = int(metaclustering.max()) + 1 if len(metaclustering) > 0 else 1
        cmap = plt.cm.get_cmap("tab20", n_meta)
        max_sz = float(node_sizes.max()) if node_sizes.max() > 0 else 1.0
        display_sizes = 80 + (node_sizes / max_sz) * 500

        fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)
        ax.set_facecolor(BG_COLOR)

        # Arêtes
        for i, j in edges:
            ax.plot(
                [layout_coords[i, 0], layout_coords[j, 0]],
                [layout_coords[i, 1], layout_coords[j, 1]],
                color=SPINE_COLOR,
                lw=1.2,
                zorder=1,
                alpha=0.7,
            )

        # Nodes
        node_colors = [
            cmap(mc_per_node[i] % 20)
            if mc_per_node[i] >= 0
            else (0.45, 0.45, 0.45, 1.0)
            for i in range(n_nodes)
        ]
        ax.scatter(
            layout_coords[:, 0],
            layout_coords[:, 1],
            s=display_sizes,
            c=node_colors,
            zorder=2,
            edgecolors=TEXT_COLOR,
            linewidths=0.5,
            alpha=0.9,
        )

        # Labels — centroïde de chaque métacluster
        mc_pts: dict = defaultdict(list)
        for i in range(n_nodes):
            if mc_per_node[i] >= 0:
                mc_pts[mc_per_node[i]].append(layout_coords[i])
        for mc_id, pts in mc_pts.items():
            centroid = np.mean(pts, axis=0)
            ax.text(
                centroid[0],
                centroid[1],
                f"MC{mc_id}",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color="white",
                zorder=3,
            )

        ax.set_title(title, fontsize=14, fontweight="bold", color=TEXT_COLOR, pad=12)
        ax.set_xlabel("MST Dim 1", fontsize=11, color=TEXT_COLOR)
        ax.set_ylabel("MST Dim 2", fontsize=11, color=TEXT_COLOR)
        ax.tick_params(colors=TEXT_COLOR)
        for spine in ax.spines.values():
            spine.set_color(SPINE_COLOR)

        # Légende métaclusters
        legend_handles = [
            plt.scatter([], [], s=80, color=cmap(i % 20), label=f"MC{i}")
            for i in range(n_meta)
        ]
        ax.legend(
            handles=legend_handles,
            loc="best",
            facecolor="#313244",
            labelcolor=TEXT_COLOR,
            edgecolor=SPINE_COLOR,
            fontsize=9,
            title="Métacluster",
            title_fontsize=9,
        )

        save_figure(fig, output_path)
        _logger.info("MST statique sauvegardé: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_mst_static: %s", exc)
        return None


# ─────────────────────────────────────────────────────────────────────────────
#  MST Plotly interactif
# ─────────────────────────────────────────────────────────────────────────────


def plot_mst_plotly(
    clusterer: Any,
    mfi_matrix: pd.DataFrame,
    metaclustering: np.ndarray,
    output_path: Optional[Path | str] = None,
    *,
    title: str = "MST FlowSOM — Vue Interactive",
) -> Optional[Any]:
    """
    MST FlowSOM interactif Plotly.

    Nodes colorés par métacluster, taille proportionnelle au nb de cellules.
    Hover: ID node, métacluster, cellules, top 2 marqueurs.

    Args:
        clusterer: FlowSOMClusterer après fit().
        mfi_matrix: DataFrame [n_metaclusters × n_markers].
        metaclustering: Assignation métacluster par cellule.
        output_path: Chemin HTML de sortie (optionnel).
        title: Titre du graphique.

    Returns:
        Figure Plotly ou None.
    """
    try:
        import plotly.graph_objects as go
        import plotly.express as px
        from scipy.spatial.distance import cdist
        from scipy.sparse.csgraph import minimum_spanning_tree
        from scipy.sparse import csr_matrix
    except ImportError:
        _logger.warning("plotly / scipy requis pour plot_mst_plotly")
        return None

    try:
        layout_coords = clusterer.get_layout_coords()
        node_sizes = clusterer.get_node_sizes()
        n_nodes = clusterer.n_nodes

        # Métacluster dominant par node
        mc_per_node = np.full(n_nodes, -1, dtype=int)
        na = getattr(clusterer, "node_assignments_", None)
        ma = getattr(clusterer, "metacluster_assignments_", None)
        if na is not None and ma is not None:
            for node_id in range(n_nodes):
                mask_node = na == node_id
                if mask_node.any():
                    mc_per_node[node_id] = int(np.bincount(ma[mask_node]).argmax())

        # Arêtes MST
        edge_x: List[Optional[float]] = []
        edge_y: List[Optional[float]] = []
        fsm = getattr(clusterer, "_fsom_model", None)
        codebook = None
        if fsm is not None:
            if hasattr(fsm, "codes"):
                codebook = np.asarray(fsm.codes, dtype=float)
            elif hasattr(fsm, "model") and hasattr(fsm.model, "codes"):
                codebook = np.asarray(fsm.model.codes, dtype=float)
        if codebook is not None:
            dist_mat = cdist(codebook, codebook, metric="euclidean")
            mst = minimum_spanning_tree(csr_matrix(dist_mat))
            coo = mst.tocoo()
            for i, j in zip(coo.row.tolist(), coo.col.tolist()):
                edge_x += [layout_coords[i, 0], layout_coords[j, 0], None]
                edge_y += [layout_coords[i, 1], layout_coords[j, 1], None]

        palette = (
            px.colors.qualitative.Set1
            + px.colors.qualitative.Pastel
            + px.colors.qualitative.Set2
        )
        n_meta = int(metaclustering.max()) + 1 if len(metaclustering) > 0 else 1
        max_sz = float(node_sizes.max()) if node_sizes.max() > 0 else 1.0

        traces: List[Any] = []
        # Arêtes en premier
        if edge_x:
            traces.append(
                go.Scatter(
                    x=edge_x,
                    y=edge_y,
                    mode="lines",
                    line=dict(color="rgba(180,180,180,0.35)", width=1.5),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

        # Nodes par métacluster
        for mc_id in range(n_meta):
            indices = [i for i in range(n_nodes) if mc_per_node[i] == mc_id]
            if not indices:
                continue
            idx = np.array(indices)
            node_sz = 10 + (node_sizes[idx] / max_sz) * 38

            mc_key = f"MC{mc_id}"
            if mc_key in mfi_matrix.index:
                top2 = mfi_matrix.loc[mc_key].nlargest(2).index.tolist()
                hover = [
                    f"Node {i}<br><b>MC{mc_id}</b><br>Cellules: {int(node_sizes[i]):,}<br>Top: {', '.join(top2)}"
                    for i in idx
                ]
            else:
                hover = [
                    f"Node {i}<br><b>MC{mc_id}</b><br>Cellules: {int(node_sizes[i]):,}"
                    for i in idx
                ]

            traces.append(
                go.Scatter(
                    x=layout_coords[idx, 0].tolist(),
                    y=layout_coords[idx, 1].tolist(),
                    mode="markers+text",
                    marker=dict(
                        size=node_sz.tolist(),
                        color=palette[mc_id % len(palette)],
                        line=dict(color="white", width=0.8),
                        opacity=0.9,
                    ),
                    text=[f"MC{mc_id}"] * len(idx),
                    textposition="middle center",
                    textfont=dict(size=8, color="white"),
                    name=f"MC{mc_id}",
                    hovertext=hover,
                    hoverinfo="text",
                )
            )

        fig = go.Figure(data=traces)
        fig.update_layout(
            title=dict(text=f"<b>{title}</b>", font=dict(size=16)),
            paper_bgcolor="#1e1e2e",
            plot_bgcolor="#1e1e2e",
            font=dict(color="#e2e8f0"),
            height=620,
            showlegend=True,
            legend=dict(bgcolor="#313244", bordercolor="#585b70"),
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
            margin=dict(l=20, r=20, t=60, b=20),
        )

        if output_path is not None:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(out), include_plotlyjs="cdn")
            _logger.info("MST Plotly sauvegardé: %s", out.name)

        return fig

    except Exception as exc:
        _logger.error("Échec plot_mst_plotly: %s", exc)
        return None


# ─────────────────────────────────────────────────────────────────────────────
#  Grille SOM Plotly (ScatterGL)
# ─────────────────────────────────────────────────────────────────────────────


def plot_som_grid_plotly(
    clustering: np.ndarray,
    metaclustering: np.ndarray,
    clusterer: Any,
    output_path: Optional[Path | str] = None,
    *,
    max_cells: int = 50_000,
    seed: int = 42,
    title: str = "Grille SOM — Cellules par Métacluster",
) -> Optional[Any]:
    """
    Grille SOM Plotly ScatterGL — chaque cellule positionnée sur son node SOM avec jitter.

    Args:
        clustering: Assignation de node par cellule (n_cells,).
        metaclustering: Assignation métacluster par cellule (n_cells,).
        clusterer: FlowSOMClusterer après fit().
        output_path: Chemin HTML de sortie (optionnel).
        max_cells: Sous-échantillonnage si > max_cells cellules.
        seed: Graine pour reproducibilité.
        title: Titre du graphique.

    Returns:
        Figure Plotly ou None.
    """
    try:
        import plotly.graph_objects as go
        import plotly.express as px
    except ImportError:
        _logger.warning("plotly requis pour plot_som_grid_plotly")
        return None

    try:
        grid_coords = clusterer.get_grid_coords()  # (n_nodes, 2)
        node_sizes = clusterer.get_node_sizes()  # (n_nodes,)

        # Sous-échantillonnage
        rng = np.random.default_rng(seed)
        n = min(max_cells, len(clustering))
        idx_sample = rng.choice(len(clustering), n, replace=False)
        cl_sub = clustering[idx_sample].astype(int)
        mc_sub = metaclustering[idx_sample].astype(int)

        # Jitter circulaire proportionnel à la taille du node
        max_sz = float(node_sizes.max()) if node_sizes.max() > 0 else 1.0
        theta = rng.uniform(0.0, 2.0 * np.pi, n)
        u = rng.uniform(0.0, 1.0, n)
        radii = 0.1 + 0.35 * np.sqrt(node_sizes[cl_sub] / max_sz)
        r = np.sqrt(u) * radii
        x_pos = grid_coords[cl_sub, 0] + r * np.cos(theta)
        y_pos = grid_coords[cl_sub, 1] + r * np.sin(theta)

        palette = (
            px.colors.qualitative.Set1
            + px.colors.qualitative.Pastel
            + px.colors.qualitative.Set2
        )
        n_meta = int(metaclustering.max()) + 1 if len(metaclustering) > 0 else 1

        traces: List[Any] = []
        for mc_id in range(n_meta):
            mask = mc_sub == mc_id
            if not mask.any():
                continue
            traces.append(
                go.Scattergl(
                    x=x_pos[mask].tolist(),
                    y=y_pos[mask].tolist(),
                    mode="markers",
                    marker=dict(
                        size=3,
                        color=palette[mc_id % len(palette)],
                        opacity=0.55,
                    ),
                    name=f"MC{mc_id}",
                    hoverinfo="skip",
                )
            )

        fig = go.Figure(data=traces)
        fig.update_layout(
            title=dict(text=f"<b>{title}</b>", font=dict(size=16)),
            paper_bgcolor="#1e1e2e",
            plot_bgcolor="#1e1e2e",
            font=dict(color="#e2e8f0"),
            height=620,
            legend=dict(bgcolor="#313244", bordercolor="#585b70"),
            xaxis=dict(
                title="SOM Grid X",
                showgrid=False,
                zeroline=False,
                color="#e2e8f0",
            ),
            yaxis=dict(
                title="SOM Grid Y",
                showgrid=False,
                zeroline=False,
                color="#e2e8f0",
            ),
            margin=dict(l=40, r=20, t=60, b=40),
        )

        if output_path is not None:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(out), include_plotlyjs="cdn")
            _logger.info("SOM Grid Plotly sauvegardé: %s", out.name)

        return fig

    except Exception as exc:
        _logger.error("Échec plot_som_grid_plotly: %s", exc)
        return None


# ─────────────────────────────────────────────────────────────────────────────
#  Optimisation FlowSOM — Visualisation multi-critères
# ─────────────────────────────────────────────────────────────────────────────


def plot_optimization_results(
    results_df: pd.DataFrame,
    best_k: int,
    stability_results: Optional[Dict[int, Any]] = None,
    w_stability: float = 0.65,
    w_silhouette: float = 0.35,
    min_stability_threshold: float = 0.75,
    output_path: Optional[Path | str] = None,
) -> Optional["_mpl_figure.Figure"]:
    """
    Visualisation des résultats d'optimisation multi-critères FlowSOM.

    Produit une figure matplotlib à 2 ou 3 panneaux :
      - Panneau 1 : Silhouette score vs k
      - Panneau 2 : ARI moyen/stabilité bootstrap (si stability_results fourni)
      - Panneau 3 (ou 2) : Score composite pondéré (w_stability × ARI + w_silhouette × Sil)

    Args:
        results_df: DataFrame avec colonnes 'k', 'silhouette', optionnellement
                    'composite_score'.
        best_k: Valeur optimale de k retenue.
        stability_results: Dict {k: {'mean_ari': float, 'std_ari': float}} (optionnel).
        w_stability: Poids de l'ARI dans le score composite (défaut 0.65).
        w_silhouette: Poids du silhouette dans le score composite (défaut 0.35).
        min_stability_threshold: Seuil ARI minimal pour déclarer un k stable.
        output_path: Chemin PNG de sauvegarde (optionnel).

    Returns:
        Figure matplotlib ou None si matplotlib absent.
    """
    if not _MPL_AVAILABLE:
        _logger.warning("matplotlib requis pour plot_optimization_results")
        return None

    has_stability = bool(stability_results and len(stability_results) > 0)
    n_plots = 3 if has_stability else 2

    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5), facecolor=BG_COLOR)
    axes = list(axes)

    ks = results_df["k"].values
    sils = results_df["silhouette"].values

    # ── Panneau 1 : Silhouette ─────────────────────────────────────────────
    ax = axes[0]
    ax.set_facecolor(BG_COLOR)
    ax.plot(
        ks, sils, "o-", color="#2196F3", linewidth=2, markersize=5, label="Silhouette"
    )
    ax.axvline(
        best_k,
        color="#F44336",
        linestyle="--",
        linewidth=2,
        alpha=0.7,
        label=f"Optimal k={best_k}",
    )
    ax.set_xlabel("Nombre de métaclusters (k)", fontsize=11, color=TEXT_COLOR)
    ax.set_ylabel("Silhouette Score", fontsize=11, color=TEXT_COLOR)
    ax.set_title(
        "Silhouette sur Codebook SOM", fontsize=12, fontweight="bold", color=TEXT_COLOR
    )
    ax.legend(fontsize=9, facecolor="#313244", labelcolor=TEXT_COLOR)
    ax.grid(True, alpha=0.3, color=SPINE_COLOR)
    ax.tick_params(colors=TEXT_COLOR)
    for sp in ax.spines.values():
        sp.set_color(SPINE_COLOR)

    # ── Panneau 2 : Stabilité ARI ──────────────────────────────────────────
    if has_stability:
        ax = axes[1]
        ax.set_facecolor(BG_COLOR)
        stab_ks = sorted(stability_results.keys())
        stab_aris = [stability_results[k]["mean_ari"] for k in stab_ks]
        stab_stds = [stability_results[k]["std_ari"] for k in stab_ks]
        ax.errorbar(
            stab_ks,
            stab_aris,
            yerr=stab_stds,
            fmt="s-",
            color="#4CAF50",
            linewidth=2,
            markersize=6,
            capsize=3,
            label="ARI moyen ± σ",
        )
        ax.axhline(
            min_stability_threshold,
            color="#FF9800",
            linestyle=":",
            linewidth=1.5,
            label=f"Seuil stabilité ({min_stability_threshold})",
        )
        ax.axvline(best_k, color="#F44336", linestyle="--", linewidth=2, alpha=0.7)
        ax.set_xlabel("Nombre de métaclusters (k)", fontsize=11, color=TEXT_COLOR)
        ax.set_ylabel("ARI moyen (stabilité)", fontsize=11, color=TEXT_COLOR)
        ax.set_title(
            "Stabilité Bootstrap (ARI pairwise)",
            fontsize=12,
            fontweight="bold",
            color=TEXT_COLOR,
        )
        ax.legend(fontsize=9, facecolor="#313244", labelcolor=TEXT_COLOR)
        ax.grid(True, alpha=0.3, color=SPINE_COLOR)
        ax.set_ylim(0, 1.05)
        ax.tick_params(colors=TEXT_COLOR)
        for sp in ax.spines.values():
            sp.set_color(SPINE_COLOR)

    # ── Panneau 3 (ou 2) : Score composite ────────────────────────────────
    ax = axes[-1]
    ax.set_facecolor(BG_COLOR)
    if "composite_score" in results_df.columns:
        valid = results_df.dropna(subset=["composite_score"])
        ax.bar(
            valid["k"],
            valid["composite_score"],
            color="#9C27B0",
            alpha=0.7,
            label="Score composite",
        )
        ax.axvline(
            best_k,
            color="#F44336",
            linestyle="--",
            linewidth=2,
            alpha=0.7,
            label=f"Optimal k={best_k}",
        )
        ax.set_xlabel("Nombre de métaclusters (k)", fontsize=11, color=TEXT_COLOR)
        ax.set_ylabel("Score composite", fontsize=11, color=TEXT_COLOR)
        ax.set_title(
            f"Score Composite (w_stab={w_stability}, w_sil={w_silhouette})",
            fontsize=12,
            fontweight="bold",
            color=TEXT_COLOR,
        )
    else:
        diffs = np.diff(sils)
        ax.plot(ks[1:], diffs, "o-", color="#FF5722", linewidth=2, markersize=4)
        ax.axvline(best_k, color="#F44336", linestyle="--", linewidth=2, alpha=0.7)
        ax.set_xlabel("k", fontsize=11, color=TEXT_COLOR)
        ax.set_ylabel("Δ Silhouette", fontsize=11, color=TEXT_COLOR)
        ax.set_title(
            "Variation Silhouette (Elbow)",
            fontsize=12,
            fontweight="bold",
            color=TEXT_COLOR,
        )
    ax.legend(fontsize=9, facecolor="#313244", labelcolor=TEXT_COLOR)
    ax.grid(True, alpha=0.3, color=SPINE_COLOR)
    ax.tick_params(colors=TEXT_COLOR)
    for sp in ax.spines.values():
        sp.set_color(SPINE_COLOR)

    plt.tight_layout()

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(out), dpi=150, bbox_inches="tight")
        _logger.info("Figure optimisation sauvegardée → %s", out)
        plt.close(fig)

    return fig


# ─────────────────────────────────────────────────────────────────────────────
#  Jitter circulaire (version non-seeded — style FlowSOM R)
# ─────────────────────────────────────────────────────────────────────────────


def circular_jitter_viz(
    n_points: int,
    cluster_ids: np.ndarray,
    node_sizes: np.ndarray,
    max_radius: float = 0.45,
    min_radius: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Jitter circulaire style FlowSOM R — vectorisé, sans graine fixe.

    Variante de ``circular_jitter`` utilisée dans les représentations SOM
    en temps réel où la reproductibilité n'est pas requise (exploration
    interactive, notebooks). Le rayon dépend proportionnellement à la
    racine carrée de la taille de chaque nœud.

    Args:
        n_points: Nombre de cellules à positionner.
        cluster_ids: Assignation de nœud par cellule (n_cells,), entiers.
        node_sizes: Taille de chaque nœud SOM (n_nodes,).
        max_radius: Rayon maximum pour le nœud le plus grand (défaut 0.45).
        min_radius: Rayon minimum pour les petits nœuds (défaut 0.10).

    Returns:
        Tuple (jitter_x, jitter_y) de dtype float32.
    """
    theta = np.random.uniform(0, 2 * np.pi, n_points)
    u = np.random.uniform(0, 1, n_points)

    max_size_val = float(node_sizes.max())
    if max_size_val <= 0.0:
        max_size_val = 1.0

    size_ratios = np.sqrt(node_sizes[cluster_ids.astype(int)] / max_size_val)
    radii = (min_radius + (max_radius - min_radius) * size_ratios).astype(np.float32)
    r = np.sqrt(u) * radii

    return (r * np.cos(theta)).astype(np.float32), (r * np.sin(theta)).astype(
        np.float32
    )


# =============================================================================
# SECTION §13 — Star Chart (MST view via fs.pl.plot_stars)
# =============================================================================


def plot_star_chart(
    fsom: Any,
    output_path: Path | str,
    *,
    title: str = "FlowSOM Star Chart (MST View)",
    dpi: int = 150,
) -> Optional["_mpl_figure.Figure"]:
    """
    Génère le Star Chart FlowSOM en vue MST via ``fs.pl.plot_stars()``.

    Le Star Chart est la visualisation officielle du package FlowSOM (R/Python) :
    chaque nœud du MST est représenté par un étoile dont les branches traduisent
    l'intensité relative de chaque marqueur dans ce nœud. La couleur de fond
    correspond au métacluster dominant du nœud.

    Args:
        fsom: Objet FlowSOM entraîné (retourné par ``fs.FlowSOM()``).
        output_path: Chemin PNG de sauvegarde.
        title: Titre de la figure.
        dpi: Résolution de sortie (défaut 150).

    Returns:
        Figure matplotlib ou None si la génération échoue.
    """
    if not _MPL_AVAILABLE:
        _logger.warning("matplotlib requis pour plot_star_chart")
        return None

    try:
        import flowsom as fs  # noqa: F401 — vérifie la disponibilité

        fig = fs.pl.plot_stars(
            fsom,
            background_values=fsom.get_cluster_data().obs.metaclustering,
            view="MST",
        )
        plt.suptitle(title, fontsize=14, fontweight="bold")
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
        plt.close("all")
        _logger.info("Star Chart sauvegardé: %s", output_path)
        return fig
    except Exception as exc:
        _logger.warning(
            "Échec plot_star_chart (fs.pl.plot_stars indisponible): %s", exc
        )
        return None


# =============================================================================
# SECTION §13b — Star Chart custom (compatible CPU + GPU)
#
# Reproduit visuellement le Star Chart FlowSOM sans dépendre de fs.pl.plot_stars:
#   • Chaque nœud est dessiné en étoile dont les branches = MFI normalisée par marqueur
#   • Taille du disque central ∝ nombre de cellules dans le nœud
#   • Couleur du disque = métacluster dominant
#   • Les arêtes MST sont tracées entre les nœuds voisins
#   • Layout Kamada-Kawai (identique à fs.FlowSOM.build_MST)
#
# Compatible GPU (utilise clusterer.get_layout_coords() + codebook.codes)
# et CPU (mêmes appels — get_layout_coords() retourne obsm["layout"] si disponible).
# =============================================================================


def plot_star_chart_custom(
    clusterer: Any,
    output_path: "Path | str",
    *,
    marker_names: "Optional[List[str]]" = None,
    title: str = "FlowSOM Star Chart (MST View — Custom)",
    dpi: int = 150,
    figsize: "Tuple[int, int]" = (14, 12),
    star_scale: float = 0.35,
    max_bubble: float = 500.0,
    min_bubble: float = 30.0,
) -> "Optional[_mpl_figure.Figure]":
    """
    Génère un Star Chart FlowSOM compatible CPU **et** GPU.

    Contrairement à ``plot_star_chart`` (qui délègue à ``fs.pl.plot_stars``),
    cette fonction construit le graphique entièrement avec matplotlib à partir
    des données exposées par ``FlowSOMClusterer`` :
      - ``clusterer.get_layout_coords()``  → coordonnées MST (n_nodes, 2)
      - ``clusterer.get_node_sizes()``     → tailles des nœuds (n_nodes,)
      - ``clusterer.metacluster_map_``     → métacluster par nœud (n_nodes,)
      - ``clusterer._fsom_model.codes``    → codebook (n_nodes, n_markers)

    Args:
        clusterer: Instance ``FlowSOMClusterer`` après ``fit()``.
        output_path: Chemin PNG de sauvegarde.
        marker_names: Noms des marqueurs (ordonnés comme le codebook).
            Si None, utilise M0, M1, …
        title: Titre de la figure.
        dpi: Résolution en sortie.
        figsize: Taille de la figure en pouces.
        star_scale: Rayon maximal d'une branche étoile (unités layout).
        max_bubble: Taille matplotlib maximale du disque central.
        min_bubble: Taille matplotlib minimale du disque central.

    Returns:
        Figure matplotlib ou None si la génération échoue.
    """
    if not _MPL_AVAILABLE:
        _logger.warning("matplotlib requis pour plot_star_chart_custom")
        return None

    try:
        # ── 1. Récupération du codebook ───────────────────────────────────────
        fsom_model = getattr(clusterer, "_fsom_model", None)
        codebook: Optional[np.ndarray] = None

        if fsom_model is not None:
            if hasattr(fsom_model, "codes"):
                codebook = np.asarray(fsom_model.codes, dtype=float)
            elif hasattr(fsom_model, "get_cluster_data"):
                # CPU fs.FlowSOM : codebook = .X de cluster_data
                try:
                    codebook = np.asarray(fsom_model.get_cluster_data().X, dtype=float)
                except Exception:
                    pass

        if codebook is None or codebook.shape[0] == 0:
            _logger.warning(
                "plot_star_chart_custom : codebook indisponible — star chart ignoré."
            )
            return None

        n_nodes, n_markers = codebook.shape

        # ── 2. Noms des marqueurs ─────────────────────────────────────────────
        if marker_names is None or len(marker_names) != n_markers:
            marker_names = [f"M{i}" for i in range(n_markers)]

        # ── 3. Coordonnées MST (layout Kamada-Kawai) ──────────────────────────
        layout = clusterer.get_layout_coords()  # (n_nodes, 2)

        # ── 4. Arêtes MST (optionnel) ─────────────────────────────────────────
        mst_edges: List[Tuple[int, int]] = []
        try:
            import igraph as ig
            from scipy.spatial.distance import cdist

            adj = cdist(codebook, codebook, metric="euclidean")
            g_full = ig.Graph.Weighted_Adjacency(adj, mode="undirected", loops=False)
            g_mst = ig.Graph.spanning_tree(g_full, weights=g_full.es["weight"])
            mst_edges = [(e.source, e.target) for e in g_mst.es]
        except Exception:
            pass  # Arêtes non disponibles — star chart sans lignes MST

        # ── 5. Tailles et métaclusters ────────────────────────────────────────
        node_sizes = clusterer.get_node_sizes()  # (n_nodes,)
        mc_map = getattr(clusterer, "metacluster_map_", None)
        if mc_map is None:
            mc_map = np.zeros(n_nodes, dtype=int)
        mc_map = np.asarray(mc_map, dtype=int)

        # Nombre de métaclusters
        n_meta = int(mc_map.max()) + 1
        cmap_mc = plt.cm.tab20 if n_meta <= 20 else plt.cm.turbo

        # ── 6. Normalisation MFI par marqueur (min-max sur tous les nœuds) ───
        col_min = codebook.min(axis=0)
        col_max = codebook.max(axis=0)
        col_range = np.where(col_max - col_min > 0, col_max - col_min, 1.0)
        codebook_norm = (codebook - col_min) / col_range  # ∈ [0, 1]

        # ── 7. Construction des angles des branches ───────────────────────────
        angles = np.linspace(0, 2 * np.pi, n_markers, endpoint=False)

        # ── 8. Taille des bulles ──────────────────────────────────────────────
        max_sz = float(node_sizes.max()) if node_sizes.max() > 0 else 1.0
        bubble_sz = min_bubble + (max_bubble - min_bubble) * (node_sizes / max_sz)

        # ── 9. Mise en page ───────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_facecolor("#1e1e2e")
        fig.patch.set_facecolor("#1e1e2e")

        x_vals = layout[:, 0]
        y_vals = layout[:, 1]

        # ── 10. Tracé des arêtes MST ──────────────────────────────────────────
        for s, t in mst_edges:
            ax.plot(
                [x_vals[s], x_vals[t]],
                [y_vals[s], y_vals[t]],
                color="rgba(200,200,200,0.3)" if False else "#666688",
                lw=1.2,
                alpha=0.45,
                zorder=1,
            )

        # ── 11. Tracé des étoiles et disques par nœud ────────────────────────
        for node_id in range(n_nodes):
            xc, yc = x_vals[node_id], y_vals[node_id]
            mc_id = int(mc_map[node_id])
            color_mc = cmap_mc(mc_id / max(n_meta - 1, 1))

            # Disque central (taille ∝ cellules)
            ax.scatter(
                xc,
                yc,
                s=bubble_sz[node_id],
                color=color_mc,
                edgecolors="white",
                linewidths=0.8,
                zorder=3,
                alpha=0.92,
            )

            # Branches en étoile — une par marqueur
            profile = codebook_norm[node_id]  # ∈ [0, 1]

            # Polygone étoile (fermeture du contour)
            star_x = [
                xc + star_scale * profile[j] * np.cos(angles[j])
                for j in range(n_markers)
            ]
            star_y = [
                yc + star_scale * profile[j] * np.sin(angles[j])
                for j in range(n_markers)
            ]
            star_x.append(star_x[0])
            star_y.append(star_y[0])

            ax.plot(star_x, star_y, color=color_mc, lw=0.7, alpha=0.65, zorder=2)
            ax.fill(star_x, star_y, color=color_mc, alpha=0.12, zorder=2)

            # Label métacluster (numéro, centré sur le disque)
            if node_sizes[node_id] > 0:
                ax.text(
                    xc,
                    yc,
                    str(mc_id),
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white",
                    fontweight="bold",
                    zorder=4,
                )

        # ── 12. Légende des marqueurs (guide des angles) ─────────────────────
        # Afficher les 8 marqueurs les plus variables (pour ne pas surcharger)
        marker_variance = codebook.var(axis=0)
        top_k = min(8, n_markers)
        top_idx = np.argsort(marker_variance)[::-1][:top_k]

        for j in top_idx:
            ax.annotate(
                marker_names[j],
                xy=(
                    layout[:, 0].mean() + (star_scale * 1.6) * np.cos(angles[j]),
                    layout[:, 1].mean() + (star_scale * 1.6) * np.sin(angles[j]),
                ),
                ha="center",
                va="center",
                fontsize=7,
                color="#cccccc",
                zorder=5,
            )

        # ── 13. Légende des métaclusters ──────────────────────────────────────
        from matplotlib.patches import Patch  # noqa: PLC0415

        legend_elements = [
            Patch(
                facecolor=cmap_mc(mc / max(n_meta - 1, 1)),
                edgecolor="white",
                label=f"MC{mc}",
            )
            for mc in range(n_meta)
        ]
        ax.legend(
            handles=legend_elements,
            loc="lower right",
            fontsize=8,
            framealpha=0.6,
            facecolor="#2a2a3e",
            edgecolor="#555",
            labelcolor="white",
            ncol=max(1, n_meta // 8),
            title="Métacluster",
            title_fontsize=8,
        )

        ax.set_title(title, fontsize=13, fontweight="bold", color="#e2e8f0", pad=10)

        # ── 14. Sauvegarde ────────────────────────────────────────────────────
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(
            str(output_path),
            dpi=dpi,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
        )
        plt.close("all")
        _logger.info("Star Chart custom sauvegardé: %s", output_path)
        return fig

    except Exception as exc:
        _logger.warning("Échec plot_star_chart_custom: %s", exc)
        return None


# =============================================================================
# SECTION §14 — SOM Grid statique matplotlib (2 panneaux: MC + Condition)
# =============================================================================


def plot_som_grid_static(
    clustering: np.ndarray,
    metaclustering: np.ndarray,
    grid_coords: np.ndarray,
    condition_labels: Optional[np.ndarray],
    xdim: int,
    ydim: int,
    output_path: Path | str,
    *,
    seed: int = 42,
    max_radius: float = 0.45,
    min_radius: float = 0.10,
    dpi: int = 150,
    title_prefix: str = "Grille FlowSOM",
) -> Optional["_mpl_figure.Figure"]:
    """
    Grille SOM statique matplotlib — 2 panneaux côte-à-côte.

    Reproduit exactement la section §14 du pipeline monolithique :
      - Panneau gauche  : cellules colorées par métacluster (tab20)
        avec labels de métacluster annotés au centre de chaque nœud.
      - Panneau droite  : cellules colorées par condition (Sain=vert / Patho=rouge)
        avec légende.

    Le jitter circulaire (style FlowSOM R) positionne les cellules autour
    de leur nœud avec un rayon proportionnel à √(taille_nœud).

    Args:
        clustering: Assignation de nœud par cellule (n_cells,), dtype int.
        metaclustering: Assignation de métacluster par nœud (n_nodes,).
        grid_coords: Coordonnées (n_nodes, 2) dans la grille SOM.
        condition_labels: Condition par cellule (n_cells,) — ``None`` désactive
            le panneau droit (figure 1×1).
        xdim: Largeur de la grille SOM.
        ydim: Hauteur de la grille SOM.
        output_path: Chemin PNG de sauvegarde.
        seed: Graine pour le jitter circulaire.
        max_radius: Rayon jitter max.
        min_radius: Rayon jitter min.
        dpi: Résolution de sortie.
        title_prefix: Préfixe pour les titres des panneaux.

    Returns:
        Figure matplotlib ou None si matplotlib absent.
    """
    if not _MPL_AVAILABLE:
        _logger.warning("matplotlib requis pour plot_som_grid_static")
        return None

    try:
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch

        n_cells = len(clustering)
        n_nodes = len(metaclustering)
        cluster_ids_int = clustering.astype(int)

        # ── Vectorisé : plus de list comprehension sur les 787k cellules ──────
        # Coordonnées de grille par cellule
        xGrid_base = grid_coords[cluster_ids_int, 0].astype(np.float32)
        yGrid_base = grid_coords[cluster_ids_int, 1].astype(np.float32)
        xGrid_shifted = xGrid_base - xGrid_base.min() + 1
        yGrid_shifted = yGrid_base - yGrid_base.min() + 1

        # Métacluster par cellule
        metaclustering_cells = metaclustering[cluster_ids_int]

        # Taille de chaque nœud (vectorisé via bincount)
        node_sizes = np.bincount(cluster_ids_int, minlength=n_nodes).astype(np.float32)

        # Jitter circulaire vectorisé (reproductible)
        np.random.seed(seed)
        jitter_x, jitter_y = circular_jitter_viz(
            n_cells,
            cluster_ids_int,
            node_sizes,
            max_radius=max_radius,
            min_radius=min_radius,
        )

        n_panels = 2 if condition_labels is not None else 1
        fig, axes = plt.subplots(
            1, n_panels, figsize=(8 * n_panels, 7), facecolor="white"
        )
        if n_panels == 1:
            axes = [axes]

        # ── Panneau 1 : Métaclusters ─────────────────────────────────────────
        ax1 = axes[0]
        n_meta = len(np.unique(metaclustering))
        cmap_mc = plt.cm.tab20 if n_meta <= 20 else plt.cm.turbo

        scatter1 = ax1.scatter(
            xGrid_shifted + jitter_x,
            yGrid_shifted + jitter_y,
            c=metaclustering_cells,
            cmap=cmap_mc,
            s=5,
            alpha=0.5,
            edgecolors="none",
        )
        for node_id in range(n_nodes):
            if node_sizes[node_id] > 0:
                x_pos = grid_coords[node_id, 0] - grid_coords[:, 0].min() + 1
                y_pos = grid_coords[node_id, 1] - grid_coords[:, 1].min() + 1
                meta_id = int(metaclustering[node_id])
                ax1.annotate(
                    str(meta_id + 1),
                    (x_pos, y_pos),
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                    color="white",
                    bbox=dict(
                        boxstyle="circle,pad=0.2",
                        facecolor=cmap_mc(meta_id / max(n_meta - 1, 1)),
                        edgecolor="white",
                        alpha=0.9,
                    ),
                )
        ax1.set_xlabel("xGrid", fontsize=12, fontweight="bold")
        ax1.set_ylabel("yGrid", fontsize=12, fontweight="bold")
        ax1.set_title(
            f"{title_prefix} — {xdim}×{ydim} nœuds\nColoré par Métacluster (style FlowSOM R)",
            fontsize=12,
            fontweight="bold",
        )
        ax1.set_xlim(0.5, xdim + 1.5)
        ax1.set_ylim(0.5, ydim + 1.5)
        ax1.set_aspect("equal")
        ax1.grid(True, alpha=0.3, linestyle="--")
        plt.colorbar(scatter1, ax=ax1, label="Métacluster")

        # ── Panneau 2 : Conditions ────────────────────────────────────────────
        if condition_labels is not None:
            ax2 = axes[1]
            # Vectorisé : numpy string ops plutôt qu'une list comprehension sur 787k cellules
            cond_arr = np.asarray(condition_labels, dtype=str)
            cond_lower = np.char.lower(cond_arr)
            condition_num = np.where(
                np.isin(cond_lower, ["sain", "healthy", "nbm", "normal"]), 0, 1
            ).astype(int)
            cmap_cond = ListedColormap(["#a6e3a1", "#f38ba8"])
            ax2.scatter(
                xGrid_shifted + jitter_x,
                yGrid_shifted + jitter_y,
                c=condition_num,
                cmap=cmap_cond,
                s=5,
                alpha=0.5,
                edgecolors="none",
            )
            ax2.set_xlabel("xGrid", fontsize=12, fontweight="bold")
            ax2.set_ylabel("yGrid", fontsize=12, fontweight="bold")
            ax2.set_title(
                f"{title_prefix} — {xdim}×{ydim} nœuds\nColoré par Condition (style FlowSOM R)",
                fontsize=12,
                fontweight="bold",
            )
            ax2.set_xlim(0.5, xdim + 1.5)
            ax2.set_ylim(0.5, ydim + 1.5)
            ax2.set_aspect("equal")
            ax2.grid(True, alpha=0.3, linestyle="--")
            legend_elements = [
                Patch(facecolor="#a6e3a1", edgecolor="white", label="Sain (NBM)"),
                Patch(facecolor="#f38ba8", edgecolor="white", label="Pathologique"),
            ]
            ax2.legend(handles=legend_elements, loc="upper right")

        plt.tight_layout()
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
        plt.close("all")
        _logger.info("Grille SOM statique sauvegardée: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_som_grid_static: %s", exc)
        return None


# =============================================================================
# SECTION §15 — Radar / Spider Chart interactif Plotly (MFI par métacluster)
# =============================================================================


def plot_metacluster_radar(
    mfi_matrix: np.ndarray,
    used_markers: List[str],
    metaclustering: np.ndarray,
    output_path: Path | str,
    *,
    n_metaclusters: Optional[int] = None,
    title: str = "Profil d'Expression par Métacluster (Radar Interactif)",
) -> Optional[Any]:
    """
    Spider / Radar chart interactif Plotly — tous les métaclusters.

    Reproduit la section §15 du pipeline monolithique. Pour chaque métacluster,
    les valeurs de MFI sont normalisées entre 0 et 1 (min-max par cluster) afin
    de rendre toutes les intensités comparables sur un axe radial commun.

    Chaque trace Plotly ``go.Scatterpolar`` correspond à un métacluster.
    Le hover affiche la MFI brute et la valeur normalisée par marqueur.

    Args:
        mfi_matrix: Matrice MFI (n_metaclusters × n_markers).
        used_markers: Liste des noms de marqueurs (len == n_markers).
        metaclustering: Assignation de métacluster par cellule, pour
            indiquer les effectifs dans la légende.
        output_path: Chemin HTML de sauvegarde.
        n_metaclusters: Nombre de métaclusters ; inféré depuis ``mfi_matrix``
            si None.
        title: Titre de la figure.

    Returns:
        ``plotly.graph_objects.Figure`` ou None si plotly absent.
    """
    try:
        import plotly.graph_objects as go
        import plotly.colors as pc
    except ImportError:
        _logger.warning("plotly requis pour plot_metacluster_radar")
        return None

    try:
        if n_metaclusters is None:
            n_metaclusters = mfi_matrix.shape[0]

        if n_metaclusters <= 10:
            _palette = pc.qualitative.Set3
        elif n_metaclusters <= 20:
            _palette = pc.qualitative.Alphabet
        else:
            _palette = [
                f"hsl({int(i * 360 / n_metaclusters)},70%,55%)"
                for i in range(n_metaclusters)
            ]

        fig = go.Figure()

        for cluster_id in range(n_metaclusters):
            values = mfi_matrix[cluster_id].copy()
            v_min, v_max = float(values.min()), float(values.max())
            values_norm = (values - v_min) / (v_max - v_min + 1e-10)

            _c = _palette[cluster_id % len(_palette)]
            _n_cells = int((metaclustering == cluster_id).sum())

            if "rgb" in str(_c):
                _fill = _c.replace(")", ",0.08)").replace("rgb", "rgba")
            else:
                _fill = "rgba(128,128,128,0.05)"

            fig.add_trace(
                go.Scatterpolar(
                    r=np.append(values_norm, values_norm[0]),
                    theta=list(used_markers) + [used_markers[0]],
                    fill="toself",
                    fillcolor=_fill,
                    opacity=0.85,
                    name=f"MC{cluster_id}  ({_n_cells:,} cells)",
                    line=dict(color=_c, width=2),
                    marker=dict(size=5),
                    customdata=np.stack(
                        [
                            np.append(
                                mfi_matrix[cluster_id], mfi_matrix[cluster_id][0]
                            ),
                            np.append(values_norm, values_norm[0]),
                        ],
                        axis=-1,
                    ),
                    hovertemplate=(
                        f"<b>MC{cluster_id}</b><br>"
                        "Marqueur: %{theta}<br>"
                        "MFI brute: %{customdata[0]:.2f}<br>"
                        "Normalisé: %{customdata[1]:.3f}<extra></extra>"
                    ),
                )
            )

        fig.update_layout(
            polar=dict(
                bgcolor="#1e1e2e",
                radialaxis=dict(
                    visible=True,
                    range=[0, 1.05],
                    tickfont=dict(size=9, color="white"),
                    gridcolor="rgba(255,255,255,0.15)",
                    linecolor="rgba(255,255,255,0.3)",
                ),
                angularaxis=dict(
                    tickfont=dict(size=10, color="white"),
                    gridcolor="rgba(255,255,255,0.15)",
                    linecolor="rgba(255,255,255,0.3)",
                ),
            ),
            showlegend=True,
            title=dict(
                text=title,
                font=dict(size=16, color="white"),
                x=0.5,
            ),
            paper_bgcolor="#1e1e2e",
            plot_bgcolor="#1e1e2e",
            font=dict(color="white"),
            legend=dict(
                bgcolor="#313244",
                bordercolor="#45475a",
                font=dict(color="white"),
            ),
            height=700,
        )

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(output_path), include_plotlyjs="cdn")
        _logger.info("Radar métaclusters sauvegardé: %s", output_path)
        return fig

    except Exception as exc:
        _logger.error("Échec plot_metacluster_radar: %s", exc)
        return None


# =============================================================================
# SECTION §15 — Analyse des clusters exclusifs (100% Patho / 100% Sain)
# =============================================================================


def compute_exclusive_clusters(
    metaclustering: np.ndarray,
    condition_labels: np.ndarray,
    n_metaclusters: int,
    *,
    patho_label: str = "Pathologique",
    sain_label: str = "Sain",
) -> Dict[str, Any]:
    """
    Identifie les métaclusters exclusivement pathologiques ou exclusivement sains.

    Reproduit la section §15 du pipeline monolithique (lignes ~6840–6931).
    Un cluster est ``exclusif`` si 100 % de ses cellules proviennent d'une
    seule condition. Ces clusters représentent des populations absentes ou
    aberrantes par rapport à la moelle normale de référence.

    Args:
        metaclustering: Assignation de métacluster par cellule (n_cells,).
        condition_labels: Label de condition par cellule (n_cells,).
        n_metaclusters: Nombre total de métaclusters.
        patho_label: Valeur de la condition pathologique.
        sain_label: Valeur de la condition saine.

    Returns:
        Dict avec les clés:
          ``patho_only``  : List[Tuple[int, int]] — (cluster_id, n_cells) 100% patho
          ``sain_only``   : List[Tuple[int, int]] — (cluster_id, n_cells) 100% sain
          ``mixed``       : List[int]             — cluster_ids partagés
          ``summary_lines``: List[str]            — lignes de rapport texte
    """
    patho_only: List[Tuple[int, int]] = []
    sain_only: List[Tuple[int, int]] = []
    mixed: List[int] = []

    cond_arr = np.asarray(condition_labels)

    for cluster_id in range(n_metaclusters):
        mask = metaclustering == cluster_id
        total = int(mask.sum())
        if total == 0:
            continue
        n_patho = int((mask & (cond_arr == patho_label)).sum())
        n_sain = int((mask & (cond_arr == sain_label)).sum())
        if n_patho == total:
            patho_only.append((cluster_id, total))
        elif n_sain == total:
            sain_only.append((cluster_id, total))
        else:
            mixed.append(cluster_id)

    summary_lines: List[str] = []
    summary_lines.append("=" * 70)
    summary_lines.append("ANALYSE DES CLUSTERS EXCLUSIFS")
    summary_lines.append("=" * 70)

    if patho_only:
        total_patho = sum(c[1] for c in patho_only)
        summary_lines.append(f"\n[!] CLUSTERS 100% PATHOLOGIQUES: {len(patho_only)}")
        summary_lines.append(f"    Métaclusters : {[c[0] for c in patho_only]}")
        summary_lines.append(f"    Total cellules: {total_patho:,}")
        summary_lines.append(
            "    → Ces clusters représentent des populations UNIQUEMENT présentes chez le patient"
        )
    else:
        summary_lines.append("\n    Aucun cluster exclusivement pathologique détecté")

    if sain_only:
        total_sain = sum(c[1] for c in sain_only)
        summary_lines.append(f"\n[!] CLUSTERS 100% SAINS: {len(sain_only)}")
        summary_lines.append(f"    Métaclusters : {[c[0] for c in sain_only]}")
        summary_lines.append(f"    Total cellules: {total_sain:,}")
        summary_lines.append(
            "    → Ces clusters représentent des populations ABSENTES chez le patient"
        )
    else:
        summary_lines.append("\n    Aucun cluster exclusivement sain détecté")

    summary_lines.append(f"\n    Clusters mixtes (partagés): {len(mixed)}")

    for line in summary_lines:
        _logger.info(line)

    return {
        "patho_only": patho_only,
        "sain_only": sain_only,
        "mixed": mixed,
        "summary_lines": summary_lines,
    }
