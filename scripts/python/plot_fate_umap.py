#!/usr/bin/env python3
"""
plot_fate_umap.py — Fig 4E (CellRank2 fate probabilities on the UMAP) and Fig 4F
(cluster-17 fate violin), from the CANONICAL absorption fates.

Source of truth: cluster17_absorption/fate_probabilities_allcells.tsv (per-cell absorption
probabilities on the object's own terminal macrostates). The two NE fates for display:
  ASCL1+  = terminal {10}
  ASCL1−  = sum of terminals {1, 7, 9}
(forward NE terminals only — NOT P(NE) vs P(non-NE), which mixes mislabeled upstream sinks;
and NOT the cr_bifurcation auto-terminal fates, which are a different, messier terminal set).

  4E: two UMAPs coloured by P(ASCL1+) and P(ASCL1−) (continuous; independent colour scales).
  4F: violin of P(ASCL1+) vs P(ASCL1−) over cluster-17 cells, annotated with the within-NE
      Wilcoxon P (recomputed here; matches cl17_fate_summary.tsv = 4.76e-23).

Usage:
  python plot_fate_umap.py \
      --h5ad velocity.h5ad \
      --fates data/forfig5/cluster17_absorption/fate_probabilities_allcells.tsv \
      --umap-key X_umap --cluster-key seurat_clusters --source 17 \
      --out figures --stem fig4 --plot-format pdf png
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402

POS_COL, NEG_COL = "#e7298a", "#1b6969"   # ASCL1+ pink (cl10 palette), ASCL1- teal (cl7 palette)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True, help="velocity object (for the UMAP coords)")
    ap.add_argument("--fates", required=True, help="fate_probabilities_allcells.tsv")
    ap.add_argument("--umap-key", default="X_umap", help="obsm key (or 'velocity_umap')")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--source", default="17", help="bridge cluster for the 4F violin")
    ap.add_argument("--ascl1-pos", nargs="+", default=["10"])
    ap.add_argument("--ascl1-neg", nargs="+", default=["1", "7", "9"])
    ap.add_argument("--cmap", default="magma")
    ap.add_argument("--within-ne-ratio", action="store_true",
                    help="4E: paint the within-NE ASCL1- share P_neg/(P_pos+P_neg) on one "
                         "diverging panel instead of two absolute panels; cells with total NE "
                         "fate < --ne-fate-min are greyed (avoids the near-zero upstream floor)")
    ap.add_argument("--ne-fate-min", type=float, default=0.1,
                    help="cells with P(ASCL1+)+P(ASCL1-) below this are greyed (low NE commitment)")
    ap.add_argument("--ratio-cmap", default="RdBu_r",
                    help="diverging cmap for the within-NE ratio (high=ASCL1-, low=ASCL1+)")
    ap.add_argument("--point-size", type=float, default=3.0)
    ap.add_argument("--terminals-panel", action="store_true",
                    help="also draw Fig 4D: the two NE terminal states (ASCL1+/- clusters) "
                         "highlighted on the UMAP")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--stem", default="fig4")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    ap.add_argument("--figwidth-mm", type=float, default=None,
                    help="author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa")
    ap.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    args = ap.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        w = (args.figwidth_mm / 25.4) * w_mult if args.figwidth_mm else default_w_in
        if args.figheight_mm:
            h = (args.figheight_mm / 25.4) * h_mult
        elif args.figwidth_mm:
            h = w * (default_h_in / default_w_in)  # width-only: preserve authored aspect (avoid skew)
        else:
            h = default_h_in
        return (w, h)

    import anndata as ad
    plt = setup_style(font_size=12)
    from scipy.stats import wilcoxon

    adata = ad.read_h5ad(args.h5ad)
    if args.umap_key not in adata.obsm:
        raise SystemExit(f"[abort] {args.umap_key} not in obsm ({list(adata.obsm)})")
    um = np.asarray(adata.obsm[args.umap_key])[:, :2]

    fp = pd.read_csv(args.fates, sep="\t", index_col=0)
    fp.index = fp.index.astype(str)
    fp.columns = fp.columns.astype(str)
    # align fates to the object's cells (barcodes)
    fp = fp.reindex(adata.obs_names.astype(str))
    n_missing = int(fp.iloc[:, 0].isna().sum())
    if n_missing:
        print(f"[warn] {n_missing}/{adata.n_obs} cells have no fate row (barcode mismatch?) "
              f"— painted grey")

    pos = [c for c in args.ascl1_pos if c in fp.columns]
    neg = [c for c in args.ascl1_neg if c in fp.columns]
    miss = [c for c in args.ascl1_pos + args.ascl1_neg if c not in fp.columns]
    if miss:
        print(f"[warn] fate columns not found, ignored: {miss} (have {list(fp.columns)})")
    if not pos or not neg:
        raise SystemExit(f"[abort] need ≥1 ASCL1+ and ASCL1- terminal column; "
                         f"pos={pos} neg={neg}")
    P_pos = fp[pos].sum(axis=1).to_numpy(dtype=float)
    P_neg = fp[neg].sum(axis=1).to_numpy(dtype=float)
    os.makedirs(args.out, exist_ok=True)

    # ---- 4D: two NE terminal states on the UMAP ----
    if args.terminals_panel:
        clD = adata.obs[args.cluster_key].astype(str).to_numpy()
        is_pos = np.isin(clD, [str(c) for c in args.ascl1_pos])
        is_neg = np.isin(clD, [str(c) for c in args.ascl1_neg])
        rest = ~(is_pos | is_neg)
        fig, ax = plt.subplots(figsize=figsize(6.0, 5.2))
        ax.scatter(um[rest, 0], um[rest, 1], s=args.point_size * 0.6, c="#dddddd",
                   lw=0, rasterized=True)
        ax.scatter(um[is_neg, 0], um[is_neg, 1], s=args.point_size, c=NEG_COL, lw=0,
                   rasterized=True, label=f"ASCL1-\nterminal {{{','.join(args.ascl1_neg)}}}")
        ax.scatter(um[is_pos, 0], um[is_pos, 1], s=args.point_size, c=POS_COL, lw=0,
                   rasterized=True, label=f"ASCL1+\nterminal {{{','.join(args.ascl1_pos)}}}")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=11); ax.set_ylabel("UMAP2", fontsize=11)
        ax.set_title("NE terminal states", fontsize=12)
        ax.legend(loc="best", fontsize=8, frameon=False, markerscale=2.5)
        for s in ax.spines.values():
            s.set_visible(False)
        savefig(fig, args.out, f"{args.stem}D_terminals_umap", args.plot_format)
        plt.close(fig)

    # ---- 4E: fate on the UMAP ----
    if args.within_ne_ratio:
        # within-NE ASCL1- share, greying cells with little NE fate (avoids the upstream
        # near-zero floor / never shows the mislabeled upstream sinks)
        denom = P_pos + P_neg
        committed = np.isfinite(denom) & (denom >= args.ne_fate_min)
        ratio = np.full(denom.shape, np.nan)
        ratio[committed] = P_neg[committed] / denom[committed]
        fig, ax = plt.subplots(figsize=figsize(6.0, 5.2))
        ax.scatter(um[~committed, 0], um[~committed, 1], s=args.point_size * 0.6, c="#dddddd",
                   lw=0, rasterized=True)                       # grey = low NE commitment
        o = np.where(committed)[0]
        sca = ax.scatter(um[o, 0], um[o, 1], s=args.point_size, c=ratio[o],
                         cmap=args.ratio_cmap, vmin=0, vmax=1, lw=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=11); ax.set_ylabel("UMAP2", fontsize=11)
        ax.set_title(f"within-NE fate (ASCL1- share)\n"
                     f"grey = NE fate < {args.ne_fate_min:g}; n coloured = {int(committed.sum()):,}",
                     fontsize=11)
        for s in ax.spines.values():
            s.set_visible(False)
        cb = fig.colorbar(sca, ax=ax, shrink=0.6, aspect=12, pad=0.02)
        cb.set_label("P(ASCL1-) / [P(ASCL1+)+P(ASCL1-)]\n(0 = ASCL1+, 1 = ASCL1-)", fontsize=8)
        cb.ax.tick_params(labelsize=8)
        savefig(fig, args.out, f"{args.stem}E_withinNE_ratio_umap", args.plot_format)
    else:
        # two absolute panels (independent colour scales)
        fig, axes = plt.subplots(1, 2, figsize=figsize(9.4, 4.4))
        for ax, vals, name in ((axes[0], P_pos, f"P(ASCL1+ {{{'+'.join(pos)}}})"),
                               (axes[1], P_neg, f"P(ASCL1- {{{'+'.join(neg)}}})")):
            finite = np.isfinite(vals)
            ax.scatter(um[~finite, 0], um[~finite, 1], s=args.point_size * 0.6, c="#dddddd",
                       lw=0, rasterized=True)
            o = np.argsort(np.where(finite, vals, -1))      # low first, high on top
            vmax = float(np.nanpercentile(vals[finite], 99)) if finite.any() else 1.0
            sca = ax.scatter(um[o, 0], um[o, 1], s=args.point_size, c=vals[o], cmap=args.cmap,
                             vmin=0, vmax=max(vmax, 1e-3), lw=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_title(name, fontsize=12)
            ax.set_xlabel("UMAP1", fontsize=11)
            for s in ax.spines.values():
                s.set_visible(False)
            cb = fig.colorbar(sca, ax=ax, shrink=0.6, aspect=12, pad=0.02)
            cb.set_label("fate probability", fontsize=9); cb.ax.tick_params(labelsize=8)
        axes[0].set_ylabel("UMAP2", fontsize=11)
        savefig(fig, args.out, f"{args.stem}E_fate_umap", args.plot_format)

    # ---- 4F: cluster-17 fate violin (within-NE) ----
    cl = adata.obs[args.cluster_key].astype(str).to_numpy()
    m = (cl == str(args.source)) & np.isfinite(P_pos) & np.isfinite(P_neg)
    pp, pn = P_pos[m], P_neg[m]
    if m.sum() == 0:
        raise SystemExit(f"[abort] no cells with {args.cluster_key}=={args.source} + fates")
    try:
        W, pval = wilcoxon(pn, pp, alternative="greater")
        pstr = f"Wilcoxon P(ASCL1- > ASCL1+) = {pval:.2e}"
    except ValueError as e:
        pstr = f"Wilcoxon n/a ({e})"

    fig, ax = plt.subplots(figsize=figsize(3.4, 4.0))
    parts = ax.violinplot([pp, pn], positions=[0, 1], showmedians=True, widths=0.8)
    for b, col in zip(parts["bodies"], (POS_COL, NEG_COL)):
        b.set_facecolor(col); b.set_alpha(0.7); b.set_edgecolor("k"); b.set_linewidth(0.4)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        if key in parts:
            parts[key].set_color("k"); parts[key].set_linewidth(0.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["ASCL1+\n{%s}" % "+".join(pos),
                                               "ASCL1-\n{%s}" % "+".join(neg)])
    ax.set_ylabel("fate probability (cluster-17 cells)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"cluster {args.source} fate (n={int(m.sum())})\n"
                 f"median ASCL1+ {np.median(pp):.3f} / ASCL1- {np.median(pn):.3f}\n{pstr}",
                 fontsize=12)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    savefig(fig, args.out, f"{args.stem}F_cluster{args.source}_fate_violin", args.plot_format)

    # ---- summary tsv ----
    pd.DataFrame([{
        "n_cluster17": int(m.sum()),
        "median_P_ASCL1pos": float(np.median(pp)),
        "median_P_ASCL1neg": float(np.median(pn)),
        "withinNE_frac_ASCL1pos": float(np.median(pp) / (np.median(pp) + np.median(pn))),
        "wilcoxon_P_neg_gt_pos": pstr,
    }]).to_csv(os.path.join(args.out, f"{args.stem}EF_fate_summary.tsv"), sep="\t", index=False)
    print(f"[ok] cluster {args.source}: median ASCL1+ {np.median(pp):.3f} / "
          f"ASCL1- {np.median(pn):.3f}  ({pstr})")
    print(f"[ok] wrote {args.stem}E_fate_umap + {args.stem}F_cluster{args.source}_fate_violin "
          f"(+ summary) to {args.out}/")


if __name__ == "__main__":
    main()
