#!/usr/bin/env python3
"""
plot_clinical_epithelial_qc.py — Supplementary Fig. S19 (Gao) / S20 (Li)

Epithelial-selection QC for the two clinical NEPC cohorts, computed *entirely from the
retained epithelial subset* (we do not possess the removed non-epithelial cells, so we
make no before/after purification claim — see supplementary methods).

What it shows, per patient, for each clinical cohort:
  Panel A  Lineage-identity dot plot: epithelial markers high, immune/endothelial/
           fibroblast markers near-absent  -> the retained cells are epithelial.
  Panel B  Disease-axis dot plot: AR/luminal, basal, neuroendocrine  -> the retained
           cells span the AD->NE axis being mapped onto the LTL331 trajectory.
  Panel C  Contamination bar: % of retained cells that score non-epithelial
           (foreign-lineage module score > epithelial module score).  This *is* the
           "residual non-epithelial <10%" claim, shown rather than asserted.

Dot color  = mean lognorm expression in the group (raw units, NOT z-scored, so genuinely
             absent markers render dark/small and the absence is honest).
Dot size   = fraction of cells in the group expressing the gene (lognorm > 0).

Scores are computed on .raw lognorm (the LTL331 object stores scaled data in .X and
lognorm in .raw; the merge preserved this).

Usage
  python plot_clinical_epithelial_qc.py \
      --h5ad harmony.annotated.h5ad \
      --out figures/supp/clinical_epithelial_qc \
      --plot-format pdf png
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _utils import setup_style
    setup_style()
except Exception:
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                                "image.interpolation": "nearest"})

# ---- marker panels ---------------------------------------------------------
LINEAGE = {
    "Epithelial":       ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "Immune":           ["PTPRC", "CD3D", "CD68", "LYZ"],
    "Endothelial":      ["PECAM1", "VWF"],
    # mesenchymal genes are EMT-expressed by tumour cells (this study's thesis),
    # so they are DISPLAYED for transparency but NOT counted as contamination.
    "Mesenchymal/EMT":  ["COL1A1", "DCN", "ACTA2"],
}
AXIS = {
    "AR / luminal":     ["AR", "KLK3", "KLK2", "FOLH1", "NKX3-1"],
    "Basal":            ["KRT5", "KRT14", "TP63"],
    "Neuroendocrine":   ["CHGA", "SYP", "NCAM1", "ASCL1", "NEUROD1", "INSM1"],
}
# Contamination = definitive NON-TUMOUR lineages only. Immune and endothelial cells
# never express EPCAM/keratins, so their markers are unambiguous. Mesenchymal genes
# are excluded (EMT tumour cells express them — counting them would contradict the
# EMT-gateway finding and manufacture false "non-epithelial" calls).
CONTAM_MARKERS = ["PTPRC", "CD3D", "CD68", "LYZ", "PECAM1", "VWF"]
CONTAM_THRESH = 1.0  # lognorm; cleanly separates true immune/endothelial cells

# disease-axis rank for ordering patients within a cohort (low -> high NE).
# Returns (NE-level, castration-level) so the x-axis runs least- to most-NE, with
# HSPC before CRPC at the same NE-level:
#   PIN(0) < adeno/PRAD/HSPC(1) < mixed-NEPC(2) < pure small-cell / NEPC(3)
def phenotype_rank(pheno: str):
    p = (pheno or "").lower()
    if "pin" in p:
        ne = 0
    elif "puresmall" in p or ("smallcell" in p and "mixed" not in p):
        ne = 3                      # Gao_NEPC_smallcell, Li pureSmallCell
    elif "mixed" in p and ("nepc" in p or "small" in p or "ne" in p):
        ne = 2                      # mixed adeno/NE biopsies
    elif "nepc" in p or "neuroendocrine" in p:
        ne = 3
    else:
        ne = 1                      # adeno / prad / hspc
    cast = 1 if "crpc" in p else 0  # HSPC before CRPC at the same NE-level
    return (ne, cast)


def genes_present(adata, genes):
    have = [g for g in genes if g in adata.raw.var_names] if adata.raw is not None \
        else [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in have]
    return have, missing


def raw_matrix(adata, genes):
    """Return dense lognorm matrix (cells x genes) for `genes` from .raw."""
    src = adata.raw if adata.raw is not None else adata
    idx = [src.var_names.get_loc(g) for g in genes]
    X = src[:, genes].X if adata.raw is None else src.X[:, idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return X


def dot_stats(adata, marker_dict, groupby, groups):
    """For each (gene, group): mean lognorm and fraction expressing.

    `groups` is the explicit, disease-axis-ordered list of group labels (plain
    strings) — passed in so ordering never depends on categorical dtype.
    """
    flat, panel_of = [], {}
    for panel, genes in marker_dict.items():
        have, missing = genes_present(adata, genes)
        if missing:
            print(f"  [warn] {panel}: not in matrix -> {missing}", file=sys.stderr)
        for g in have:
            flat.append(g)
            panel_of[g] = panel
    X = raw_matrix(adata, flat)
    obs_grp = np.asarray(adata.obs[groupby].values).astype(str)
    means = np.zeros((len(flat), len(groups)))
    fracs = np.zeros((len(flat), len(groups)))
    for j, grp in enumerate(groups):
        mask = obs_grp == grp
        sub = X[mask]
        means[:, j] = sub.mean(0)
        fracs[:, j] = (sub > 0).mean(0)
    return flat, panel_of, groups, means, fracs


BLOCK_ABBR = {
    "Epithelial": "Epi", "Immune": "Immune", "Endothelial": "Endo",
    "Mesenchymal/EMT": "Mes/\nEMT", "AR / luminal": "AR/lum",
    "Basal": "Basal", "Neuroendocrine": "NE",
}


def draw_dotplot(ax, genes, panel_of, groups, means, fracs, title, show_xlabels=True):
    # color scaled per-gene across groups is misleading for absent markers;
    # use raw mean lognorm on a shared scale so absence reads as dark+small.
    vmax = np.percentile(means, 99) if means.size else 1.0
    vmax = max(vmax, 1e-6)
    ng, ngrp = len(genes), len(groups)
    for i in range(ng):
        for j in range(ngrp):
            size = 10 + 240 * fracs[i, j]
            ax.scatter(j, ng - 1 - i, s=size, c=[means[i, j]], cmap="Reds",
                       vmin=0, vmax=vmax, edgecolors="0.55", linewidths=0.3,
                       zorder=3)
    ax.set_xticks(range(ngrp))
    if show_xlabels:
        ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=7)
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes[::-1], fontsize=10, style="italic")
    ax.set_xlim(-0.6, ngrp - 0.4)
    ax.set_ylim(-0.6, ng - 0.4)
    # panel labels in a rotated LEFT gutter (clear of gene names AND the colorbar),
    # centred on each marker block; light boundary lines between blocks.
    import matplotlib.transforms as mtransforms
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    blocks, prev = {}, None
    for i, g in enumerate(genes):
        blocks.setdefault(panel_of[g], []).append(ng - 1 - i)  # data-y of each row
    for name, ys in blocks.items():
        ax.text(-0.30, float(np.mean(ys)), BLOCK_ABBR.get(name, name),
                transform=trans, rotation=90, va="center", ha="center",
                fontsize=6.5, fontweight="bold", color="0.25", linespacing=0.9)
    boundaries = []
    for i, g in enumerate(genes):
        if prev is not None and panel_of[g] != prev:
            boundaries.append(ng - i - 0.5)
        prev = panel_of[g]
    for b in boundaries:
        ax.axhline(b, color="0.85", lw=0.6, zorder=1)
    ax.set_title(title, fontsize=9, fontweight="bold", pad=8)
    ax.tick_params(length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap="Reds",
                               norm=plt.Normalize(vmin=0, vmax=vmax))
    cb = plt.colorbar(sm, ax=ax, fraction=0.022, pad=0.03)
    cb.set_label("mean lognorm", fontsize=6)
    cb.ax.tick_params(labelsize=6)


def contamination(adata, groupby, groups):
    """% of retained cells expressing a definitive immune/endothelial marker.

    A direct, unarguable non-tumour call: immune and endothelial cells never express
    EPCAM/keratins, so detecting their markers (lognorm >= CONTAM_THRESH) flags a
    genuine non-epithelial cell. Mesenchymal genes are deliberately excluded — EMT
    tumour cells express them, so they are not evidence of contamination.
    Rows are returned in the given `groups` (disease-axis) order.
    """
    markers, missing = genes_present(adata, CONTAM_MARKERS)
    if missing:
        print(f"  [warn] contamination markers not in matrix -> {missing}",
              file=sys.stderr)
    X = raw_matrix(adata, markers)                 # cells x markers, lognorm
    is_contam = (X >= CONTAM_THRESH).any(axis=1)
    grp = np.asarray(adata.obs[groupby].values).astype(str)
    rows = []
    for g in groups:
        m = grp == g
        rows.append({"group": g, "pct_nonepithelial": 100 * is_contam[m].mean(),
                     "n_cells": int(m.sum())})
    return pd.DataFrame(rows)


def draw_contam(ax, contam_df):
    y = np.arange(len(contam_df))
    vals = contam_df["pct_nonepithelial"].values
    ax.barh(y, vals, color=np.where(vals <= 10, "#4C78A8", "#E45756"),
            edgecolor="0.3", linewidth=0.4)
    ax.axvline(10, color="0.4", ls="--", lw=0.8)
    ax.set_yticks(y)
    # patient ID only (first line) — phenotype is already on the dot-plot x-axis
    short = [str(g).split("\n")[0] for g in contam_df["group"]]
    ax.set_yticklabels(short, fontsize=10)
    ax.set_xlabel("% cells expressing\nimmune / endothelial markers", fontsize=7)
    ax.set_title("Non-tumour contamination", fontsize=9, fontweight="bold")
    for yi, v in zip(y, vals):
        ax.text(v + 0.3, yi, f"{v:.1f}%", va="center", fontsize=6.5)
    ax.set_xlim(0, max(12, vals.max() * 1.25))
    ax.invert_yaxis()  # first patient on top, matching the dot-plot's left-to-right
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)


def build_cohort_figure(adata, cohort, outbase, fmts):
    sub = adata[adata.obs["dataset"] == cohort].copy()
    if sub.n_obs == 0:
        print(f"  [skip] no cells for dataset={cohort}", file=sys.stderr)
        return None
    # order patients along the disease axis
    pat = sub.obs[["orig.ident", "patient_phenotype"]].drop_duplicates().copy()
    # apply over plain strings: mapping a CATEGORICAL column to tuples makes pandas
    # build a MultiIndex of categories and crash ("isna not defined for MultiIndex").
    ranks = [phenotype_rank(str(p)) for p in pat["patient_phenotype"]]
    pat["ne"] = [r[0] for r in ranks]
    pat["cast"] = [r[1] for r in ranks]
    pat["orig.ident"] = pat["orig.ident"].astype(str)
    pat = pat.sort_values(["ne", "cast", "orig.ident"])
    order = list(pat["orig.ident"])
    label = {r["orig.ident"]: f"{r['orig.ident']}\n{r['patient_phenotype']}"
             for _, r in pat.iterrows()}
    group_order = [label[o] for o in order]     # disease-axis-ordered display labels
    # plain-string group column (no categorical dtype -> no index/isna pitfalls)
    sub.obs["_grplab"] = sub.obs["orig.ident"].astype(str).map(label).astype(str)

    g1, po1, grp1, m1, f1 = dot_stats(sub, LINEAGE, "_grplab", group_order)
    g2, po2, grp2, m2, f2 = dot_stats(sub, AXIS, "_grplab", group_order)
    contam = contamination(sub, "_grplab", group_order)

    fig = plt.figure(figsize=(max(7, 1.1 * len(order) + 4), 9))
    gs = GridSpec(2, 2, height_ratios=[len(g1), len(g2)],
                  width_ratios=[3, 1.15], hspace=0.18, wspace=0.55, figure=fig)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[:, 1])
    # top plot hides x-labels (identical to bottom); bottom carries them once
    draw_dotplot(axA, g1, po1, grp1, m1, f1, f"{cohort}  —  lineage identity",
                 show_xlabels=False)
    draw_dotplot(axB, g2, po2, grp2, m2, f2, f"{cohort}  —  disease axis")
    draw_contam(axC, contam)
    # legend for dot size
    for frac, xo in [(0.25, 0.0), (0.5, 1.0), (1.0, 2.2)]:
        axB.scatter([], [], s=10 + 240 * frac, c="0.5", edgecolors="0.4",
                    label=f"{int(frac*100)}%")
    axB.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.02, -0.15),
               fontsize=6, title_fontsize=6, frameon=False, ncol=3)
    fig.suptitle(f"Epithelial-selection QC — {cohort} cohort",
                 fontsize=11, fontweight="bold")

    os.makedirs(os.path.dirname(os.path.abspath(outbase)), exist_ok=True)
    for fmt in fmts:
        path = f"{outbase}_{cohort}.{fmt}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  [ok] wrote {path}")
    plt.close(fig)
    return contam.assign(cohort=cohort)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", default="figures/supp/clinical_epithelial_qc")
    ap.add_argument("--cohorts", nargs="+", default=["Gao", "Li"])
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    print("[version] plot_clinical_epithelial_qc rev6 (no-categorical + flatten)")
    print(f"[load] {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    if adata.raw is None:
        print("[warn] .raw is None; scoring from .X (verify it is lognorm!)",
              file=sys.stderr)

    # --- defensive: flatten any MultiIndex on obs/var/raw.var --------------
    # (a MultiIndex trips pandas' isna checks: "isna is not defined for MultiIndex")
    def _flat(ix):
        return pd.Index(["|".join(map(str, t)) for t in ix])
    print(f"[diag] obs.index={type(adata.obs.index).__name__}, "
          f"var.index={type(adata.var.index).__name__}, "
          f"raw.var.index="
          f"{type(adata.raw.var.index).__name__ if adata.raw is not None else 'NA'}")
    if isinstance(adata.obs.index, pd.MultiIndex):
        adata.obs.index = _flat(adata.obs.index)
    if isinstance(adata.var.index, pd.MultiIndex):
        adata.var.index = _flat(adata.var.index)
    if adata.raw is not None and isinstance(adata.raw.var.index, pd.MultiIndex):
        rawad = adata.raw.to_adata()
        rawad.var.index = _flat(rawad.var.index)
        adata.raw = rawad
    # also flatten any MultiIndex column axis on obs
    if isinstance(adata.obs.columns, pd.MultiIndex):
        adata.obs.columns = _flat(adata.obs.columns)
    all_contam = []
    for cohort in args.cohorts:
        print(f"[cohort] {cohort}")
        c = build_cohort_figure(adata, cohort, args.out, args.plot_format)
        if c is not None:
            all_contam.append(c)
    if all_contam:
        tsv = f"{args.out}_contamination.tsv"
        pd.concat(all_contam).to_csv(tsv, sep="\t", index=False)
        print(f"[ok] wrote {tsv}")
        # one-line summary for the methods sentence
        allc = pd.concat(all_contam)
        print(f"[summary] max % non-epithelial across patients = "
              f"{allc['pct_nonepithelial'].max():.1f}%  "
              f"(n={len(allc)} patients)")


if __name__ == "__main__":
    main()
