"""
analyze_results.py
Genere toutes les statistiques et heatmaps depuis results_grid_search.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent
CSV_PATH     = PROJECT_ROOT / "results_grid_search.csv"
OUT_DIR      = PROJECT_ROOT / "analysis_results"
OUT_DIR.mkdir(exist_ok=True)

plt.rcParams.update({"figure.dpi": 150, "font.size": 9})

# --- Chargement --------------------------------------------------------------
df = pd.read_csv(CSV_PATH)
# Garder seulement les combos complets (tous les patients presents)
counts = df.groupby(["ratio_nbm","grid_k","max_normal_pct","min_patho_pct"])["patient_id"].count()
complete = counts[counts == counts.max()].index
df = df.set_index(["ratio_nbm","grid_k","max_normal_pct","min_patho_pct"]).loc[complete].reset_index()

n_patients  = df["patient_id"].nunique()
n_combos    = df.groupby(["ratio_nbm","grid_k","max_normal_pct","min_patho_pct"]).ngroups
print(f"CSV charge : {len(df)} lignes | {n_patients} patients | {n_combos} combos complets")

def col_k(k):
    d = int(round(k**0.5))
    return f"{d}×{d}\n({k})"

# =============================================================================
# 1. STATISTIQUES GLOBALES
# =============================================================================
print("\n" + "="*60)
print("STATISTIQUES GLOBALES")
print("="*60)

best_combo = (df.groupby(["ratio_nbm","grid_k","max_normal_pct","min_patho_pct"])["abs_error"]
              .mean().reset_index().rename(columns={"abs_error":"mae"}).sort_values("mae"))

print("\nTOP 15 configurations (MAE moyenne sur tous patients) :")
print(best_combo.head(15).to_string(index=False))

best = best_combo.iloc[0]
print(f"\n* MEILLEURE CONFIG :")
print(f"   ratio_nbm      = {best['ratio_nbm']}")
print(f"   grid_k         = {best['grid_k']} ({int(round(best['grid_k']**0.5))}×{int(round(best['grid_k']**0.5))})")
print(f"   max_normal_pct = {best['max_normal_pct']}%")
print(f"   min_patho_pct  = {best['min_patho_pct']}%")
print(f"   MAE            = {best['mae']:.3f}%")

# Correlation MRD predit vs vrai sur la meilleure config
best_df = df[
    (df["ratio_nbm"]      == best["ratio_nbm"]) &
    (df["grid_k"]         == best["grid_k"]) &
    (df["max_normal_pct"] == best["max_normal_pct"]) &
    (df["min_patho_pct"]  == best["min_patho_pct"])
]
r, p = stats.pearsonr(best_df["mrd_true"], best_df["mrd_predicted"])
rmse = np.sqrt(((best_df["mrd_true"] - best_df["mrd_predicted"])**2).mean())
print(f"\n   Pearson r = {r:.4f} (p={p:.2e})")
print(f"   RMSE      = {rmse:.3f}%")
print(f"   MAE       = {best_df['abs_error'].mean():.3f}%")
print(f"   Mediane   = {best_df['abs_error'].median():.3f}%")

# =============================================================================
# 2. HEATMAP 1 — ratio × k  (moyenne sur max_normal_pct et min_patho_pct)
# =============================================================================
def make_pivot(data, row, col, val="abs_error"):
    return (data.groupby([row, col])[val].mean().reset_index()
            .pivot(index=row, columns=col, values=val)
            .sort_index(ascending=False))

def heatmap(pivot, title, path, xlabel="", ylabel="", fmt=".1f", figsize=None):
    n_c, n_r = len(pivot.columns), len(pivot.index)
    fw = figsize[0] if figsize else max(8, n_c * 1.1)
    fh = figsize[1] if figsize else max(5, n_r * 0.85)
    fig, ax = plt.subplots(figsize=(fw, fh))
    mask = pivot.isna()
    sns.heatmap(pivot, annot=True, fmt=fmt, annot_kws={"size": max(7, min(10, 80//(n_c*n_r)+1))},
                cmap="RdYlGn_r", linewidths=0.4, ax=ax, mask=mask,
                cbar_kws={"label": "Erreur Absolue Moyenne (%)"})
    # Encadrer le minimum
    vals = pivot.values.astype(float)
    mn = np.nanmin(vals)
    for r in range(n_r):
        for c in range(n_c):
            if not np.isnan(vals[r, c]) and abs(vals[r, c] - mn) < 1e-6:
                ax.add_patch(plt.Rectangle((c, r), 1, 1, fill=False, edgecolor="blue", lw=2.5))
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel(xlabel, fontsize=10); ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(axis="x", labelsize=8, rotation=45)
    ax.tick_params(axis="y", labelsize=9, rotation=0)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path.name}")

# Heatmap ratio × k
p1 = make_pivot(df, "ratio_nbm", "grid_k")
p1.columns = [col_k(c) for c in p1.columns]
p1.index   = [f"×{r}" for r in p1.index]
heatmap(p1, f"ratio × grille SOM — MAE moyenne ({n_patients} patients, tous params MRD)",
        OUT_DIR / "h1_ratio_vs_k.png", "Taille grille k", "Ratio NBM/Patient")

# =============================================================================
# 3. HEATMAP 2 — max_normal_pct × min_patho_pct  (moyenne sur ratio et k)
# =============================================================================
p2 = make_pivot(df, "max_normal_pct", "min_patho_pct")
p2.columns = [f"minP={c:.0f}%" for c in p2.columns]
p2.index   = [f"maxN={r:.1f}%" for r in p2.index]
heatmap(p2, f"Params MRD JF — MAE moyenne (tous ratios et grilles, {n_patients} patients)",
        OUT_DIR / "h2_mrd_params.png", "min_patho_cells_pct", "max_normal_marrow_pct")

# =============================================================================
# 4. HEATMAP 3 — ratio × k  avec meilleurs params MRD fixes
# =============================================================================
best_mrd = (df.groupby(["max_normal_pct","min_patho_pct"])["abs_error"]
            .mean().idxmin())
sub_bmrd = df[(df["max_normal_pct"] == best_mrd[0]) & (df["min_patho_pct"] == best_mrd[1])]
p3 = make_pivot(sub_bmrd, "ratio_nbm", "grid_k")
p3.columns = [col_k(c) for c in p3.columns]
p3.index   = [f"×{r}" for r in p3.index]
heatmap(p3, f"ratio × grille @ maxN={best_mrd[0]:.1f}% minP={best_mrd[1]:.0f}% (meilleurs params MRD)",
        OUT_DIR / "h3_ratio_k_bestMRD.png", "Taille grille k", "Ratio NBM/Patient")

# =============================================================================
# 5. HEATMAP 4 — max_normal × min_patho  avec meilleur FlowSOM fixe
# =============================================================================
best_fs = (df.groupby(["ratio_nbm","grid_k"])["abs_error"].mean().idxmin())
sub_bfs = df[(df["ratio_nbm"] == best_fs[0]) & (df["grid_k"] == best_fs[1])]
p4 = make_pivot(sub_bfs, "max_normal_pct", "min_patho_pct")
p4.columns = [f"minP={c:.0f}%" for c in p4.columns]
p4.index   = [f"maxN={r:.1f}%" for r in p4.index]
heatmap(p4, f"Params MRD JF @ ratio={best_fs[0]} k={best_fs[1]} (meilleur FlowSOM)",
        OUT_DIR / "h4_mrd_bestFlowSOM.png", "min_patho_cells_pct", "max_normal_marrow_pct")

# =============================================================================
# 6. SCATTER — MRD predit vs MRD vrai (meilleure config globale)
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (label, sub) in zip(axes, [
    (f"Meilleure config\nratio={best['ratio_nbm']} k={best['grid_k']} maxN={best['max_normal_pct']} minP={best['min_patho_pct']}", best_df),
    (f"Tous combos (moyenne par patient)", df.groupby("patient_id")[["mrd_true","mrd_predicted","abs_error"]].mean().reset_index()),
]):
    x = sub["mrd_true"] if "mrd_true" in sub.columns else sub["mrd_true"]
    y = sub["mrd_predicted"] if "mrd_predicted" in sub.columns else sub["mrd_predicted"]
    valid = x.notna() & y.notna()
    x, y = x[valid], y[valid]

    ax.scatter(x, y, alpha=0.75, s=55, edgecolors="k", linewidths=0.4, zorder=3)
    lim = max(x.max(), y.max()) * 1.08
    ax.plot([0, lim], [0, lim], "r--", lw=1.5, label="Ideal (y=x)", zorder=2)

    # Regression lineaire
    if len(x) >= 3:
        slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
        xfit = np.linspace(0, lim, 100)
        ax.plot(xfit, slope*xfit + intercept, "b-", lw=1.5,
                label=f"Regression (r={r_val:.3f})", zorder=2)

    # Labels patients
    ids = sub["patient_id"] if "patient_id" in sub.columns else sub.index
    for xi, yi, pid in zip(x, y, ids):
        ax.annotate(str(pid).replace("BLAST110_","").replace("_P1",""),
                    (xi, yi), fontsize=6, alpha=0.7,
                    xytext=(3, 3), textcoords="offset points")

    mae_val = sub["abs_error"][valid].mean() if "abs_error" in sub.columns else np.abs(x - y).mean()
    r_v, _ = stats.pearsonr(x, y)
    rmse_v = np.sqrt(((x - y)**2).mean())
    ax.set_title(f"{label}\nMAE={mae_val:.2f}%  RMSE={rmse_v:.2f}%  r={r_v:.3f}", fontsize=9)
    ax.set_xlabel("MRD Vrai (%)", fontsize=10)
    ax.set_ylabel("MRD Predit (%)", fontsize=10)
    ax.set_xlim(-2, lim); ax.set_ylim(-2, lim)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.suptitle(f"MRD Predit vs Ground Truth — {n_patients} patients AML_Dx", fontsize=12)
plt.tight_layout()
fig.savefig(OUT_DIR / "scatter_predicted_vs_true.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: scatter_predicted_vs_true.png")

# =============================================================================
# 7. HEATMAPS PAR PATIENT (ratio × k, meilleurs params MRD)
# =============================================================================
print("\nHeatmaps par patient...")
for pid in sorted(df["patient_id"].unique()):
    sub_p = sub_bmrd[sub_bmrd["patient_id"] == pid]
    if sub_p["abs_error"].isna().all():
        continue
    pp = make_pivot(sub_p, "ratio_nbm", "grid_k")
    pp.columns = [col_k(c) for c in pp.columns]
    pp.index   = [f"×{r}" for r in pp.index]
    true_mrd = df.loc[df["patient_id"] == pid, "mrd_true"].iloc[0]
    heatmap(pp,
            f"{pid}  (True MRD = {true_mrd:.1f}%)\n@ maxN={best_mrd[0]:.1f}% minP={best_mrd[1]:.0f}%",
            OUT_DIR / f"patient_{pid}.png",
            "Taille grille k", "Ratio NBM/Patient",
            figsize=(9, 5))

# =============================================================================
# 8. RÉSUMÉ PAR PATIENT — tableau MAE meilleure config
# =============================================================================
summary_per_patient = (best_df[["patient_id","mrd_true","mrd_predicted","abs_error"]]
                       .sort_values("abs_error"))
print("\n" + "="*60)
print("RÉSUMÉ PAR PATIENT (meilleure config globale)")
print("="*60)
print(summary_per_patient.to_string(index=False))

# =============================================================================
# 9. FIGURE RÉSUMÉ — 4 panels
# =============================================================================
fig = plt.figure(figsize=(18, 12))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

# Panel A : ratio × k
ax_a = fig.add_subplot(gs[0, 0])
pivot_a = make_pivot(df, "ratio_nbm", "grid_k")
pivot_a.columns = [col_k(c) for c in pivot_a.columns]
pivot_a.index   = [f"×{r}" for r in pivot_a.index]
sns.heatmap(pivot_a, annot=True, fmt=".1f", annot_kws={"size": 8},
            cmap="RdYlGn_r", linewidths=0.3, ax=ax_a,
            cbar_kws={"label": "MAE (%)", "shrink": 0.8})
ax_a.set_title("A — ratio × grille SOM\n(tous params MRD)", fontsize=10)
ax_a.set_xlabel("Grille k"); ax_a.set_ylabel("Ratio NBM")
ax_a.tick_params(labelsize=7)

# Panel B : params MRD
ax_b = fig.add_subplot(gs[0, 1])
sns.heatmap(p2, annot=True, fmt=".1f", annot_kws={"size": 8},
            cmap="RdYlGn_r", linewidths=0.3, ax=ax_b,
            cbar_kws={"label": "MAE (%)", "shrink": 0.8})
ax_b.set_title("B — Params MRD JF\n(tous ratios et grilles)", fontsize=10)
ax_b.set_xlabel("min_patho_pct"); ax_b.set_ylabel("max_normal_pct")
ax_b.tick_params(labelsize=7)

# Panel C : scatter meilleure config
ax_c = fig.add_subplot(gs[1, 0])
ax_c.scatter(best_df["mrd_true"], best_df["mrd_predicted"],
             alpha=0.75, s=55, edgecolors="k", linewidths=0.4, zorder=3)
lim_c = max(best_df["mrd_true"].max(), best_df["mrd_predicted"].max()) * 1.08
ax_c.plot([0, lim_c], [0, lim_c], "r--", lw=1.5, label="y=x")
slope, intercept, r_val, _, _ = stats.linregress(best_df["mrd_true"], best_df["mrd_predicted"])
xf = np.linspace(0, lim_c, 100)
ax_c.plot(xf, slope*xf+intercept, "b-", lw=1.5, label=f"r={r_val:.3f}")
ax_c.set_title(f"C — Scatter meilleure config\nMAE={best['mae']:.2f}%  r={r_val:.3f}", fontsize=10)
ax_c.set_xlabel("MRD Vrai (%)"); ax_c.set_ylabel("MRD Predit (%)")
ax_c.legend(fontsize=8); ax_c.grid(alpha=0.3)

# Panel D : barplot MAE par patient (meilleure config)
ax_d = fig.add_subplot(gs[1, 1])
sp = summary_per_patient.sort_values("abs_error")
colors = ["#2ecc71" if e < 5 else "#f39c12" if e < 15 else "#e74c3c"
          for e in sp["abs_error"]]
bars = ax_d.barh(sp["patient_id"].str.replace("BLAST110_","").str.replace("_P1",""),
                  sp["abs_error"], color=colors, edgecolor="k", linewidth=0.4)
ax_d.axvline(sp["abs_error"].mean(), color="navy", linestyle="--", lw=1.5,
             label=f"MAE moy = {sp['abs_error'].mean():.1f}%")
ax_d.set_title("D — Erreur absolue par patient\n(meilleure config)", fontsize=10)
ax_d.set_xlabel("Erreur Absolue (%)"); ax_d.legend(fontsize=8)
ax_d.tick_params(axis="y", labelsize=7)

plt.suptitle(
    f"FlowSOM Grid Search — Resume ({n_patients} patients AML_Dx, {n_combos} combos)\n"
    f"Meilleure config : ratio={best['ratio_nbm']} | k={best['grid_k']} "
    f"| maxN={best['max_normal_pct']}% | minP={best['min_patho_pct']}% | MAE={best['mae']:.2f}%",
    fontsize=11, y=1.01
)
fig.savefig(OUT_DIR / "summary_4panels.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: summary_4panels.png")

print(f"\nTous les fichiers dans : {OUT_DIR}")
