#!/usr/bin/env python3
"""
plot_cl17_basal_vs_emt.py — Fig 2 panel: cluster 17 is basal-LIKE, not CLASSICAL basal addressed in
Karthaus et al. (2020) or Ireland et al. (2025)

Replaces the redundant luminal/basal cytokeratin violin once the cytokeratins live in the Fig 2C
dot plot. Single cytokeratins cannot settle "classical basal": it is a PROGRAM question. Instead, this
scores each cell for a coherent classical-basal program and an EMT/mesenchymal program and shows:

  A) basal-score x EMT-score scatter: all cells (grey density) + cluster centroids (labelled) +
     cl17 cells highlighted. cl17 sits EMT-high / basal-intermediate, away from a coherent basal.
  B) per-cluster violins (canonical adeno->NE order) for, stacked:
     classical-basal score, EMT score, TP63 expression, CAV1 expression. cl17 bar highlighted.

Honest framing: cl17 = partial basal features (modest TP63) inside a dominant EMT program —
basal-like transition state, NOT a classical KRT5/KRT14/TP63-high basal cell.

Scores via sc.tl.score_genes on .raw lognorm (LTL331 .X is scaled — see [[ltl331_h5ad_schema]]).
Outputs (in --out): fig2_basal_vs_emt.{pdf,png}, fig2_basal_emt_percluster.tsv,
fig2_basal_emt_cl17_stats.tsv.

Space-limited layout (--panels): keep only the decisive scatter in main Fig 2F and move the
per-cluster violins (basal/EMT/TP63/CAV1) to Supp S9 — run twice:
  ... --panels scatter --stem fig2F          # main 2F: the basal×EMT program-space scatter
  ... --panels scatter+tp63 --stem fig2F      # main 2F: scatter + the TP63 row
  ... --panels violins --stem figS9_basal     # Supp S9: the per-cluster violin rows
Default --panels both = the combined figure.

Usage
  python plot_cl17_basal_vs_emt.py --h5ad ltl331_base.h5ad \\
      --signatures data/cl17_basal_emt_signatures.tsv --focus-cluster 17 \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --out figures/fig2_basal_emt --panels both --plot-format pdf png
"""
import argparse
import json
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig, numkey  # noqa: E402


def load_sigs(path):
    sigs = OrderedDict()
    for line in open(path):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        sigs.setdefault(parts[0].strip(), []).append(parts[1].strip())
    return sigs


def gene_vec(adata, gene):
    src = adata.raw if adata.raw is not None else adata
    var = list(map(str, src.var_names))
    if gene not in var:
        return None
    col = src.X[:, var.index(gene)]
    return (col.toarray().ravel() if ssp.issparse(col) else np.asarray(col).ravel())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--signatures", default="data/cl17_basal_emt_signatures.tsv")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--focus-cluster", default="17")
    ap.add_argument("--callout-genes", nargs="*", default=["TP63", "CAV1"])
    ap.add_argument("--cluster-order", default=None)
    ap.add_argument("--cluster-colors", default=None)
    ap.add_argument("--exclude-clusters", nargs="*", default=["18"])
    ap.add_argument("--out", default="figures/fig2_basal_emt")
    ap.add_argument("--stem", default="fig2F",
                    help="output filename stem (figure + the 2 TSVs)")
    ap.add_argument("--panels", choices=["both", "scatter", "violins", "scatter+tp63"],
                    default="both",
                    help="which component(s) to render: 'scatter' = main 2F decisive panel; "
                         "'scatter+tp63' = scatter + ONLY the TP63 violin row "
                         "for a slightly larger main 2F; 'violins' = all per-cluster "
                         "basal/EMT/TP63/CAV1 rows for Supp S9; 'both' = combined (default). "
                         "The 2 TSVs are written in every mode.")
    ap.add_argument("--figwidth-mm", type=float, default=None,
                    help="author the panel at this FINAL width in mm (the assemble_figure --qa "
                         "onpage_mm target) so the assembler scales it ~1:1 and fonts stay readable; "
                         "the saved PDF is wider than this by any outside legend — tune with --qa")
    ap.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    ap.add_argument("--scatter-swap-axes", action="store_true",
                    help="put the 2nd program (EMT) on x and the 1st (basal) on y. EMT spans the "
                         "larger range, so on the wide x-axis the cloud fills a landscape cell "
                         "better; basal (small range) fits the short y-axis.")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    def figsize(default_w_in, default_h_in):
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh:
            return (fw / 25.4, fh / 25.4)
        if fw:   # width only -> preserve authored aspect ratio (avoids skew + empty margins)
            w = fw / 25.4
            return (w, default_h_in * w / default_w_in)
        if fh:   # height only -> preserve authored aspect ratio
            h = fh / 25.4
            return (default_w_in * h / default_h_in, h)
        return (default_w_in, default_h_in)

    adata = sc.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    if args.exclude_clusters:
        adata = adata[~adata.obs[args.cluster_key].isin(
            [str(c) for c in args.exclude_clusters])].copy()
    cats = sorted(adata.obs[args.cluster_key].unique(), key=numkey)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get("order", [])]
        order = [c for c in co if c in cats] + [c for c in cats if c not in co]
    else:
        order = cats
    colors = {}
    if args.cluster_colors:
        colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()}

    # ---- score the two programs ----
    sigs = load_sigs(args.signatures)
    use_raw = adata.raw is not None
    src = set(map(str, (adata.raw.var_names if use_raw else adata.var_names)))
    score_cols = {}
    for prog, genes in sigs.items():
        present = [g for g in genes if g in src]
        missing = [g for g in genes if g not in src]
        if missing:
            print(f"[warn] {prog}: dropped {len(missing)} absent symbol(s): {missing}")
        name = f"score_{prog}"
        sc.tl.score_genes(adata, present, score_name=name, use_raw=use_raw)
        score_cols[prog] = name
    progs = list(sigs.keys())
    basal_c, emt_c = score_cols[progs[0]], score_cols[progs[1]]

    cl = adata.obs[args.cluster_key].values
    df = adata.obs[[basal_c, emt_c]].copy()
    for g in args.callout_genes:
        v = gene_vec(adata, g)
        if v is not None:
            df[g] = v
        else:
            print(f"[warn] callout gene {g} absent — skipped")
    df["cluster"] = cl

    # ---- per-cluster means + cl17 enrichment stats ----
    percluster = df.groupby("cluster").mean(numeric_only=True).reindex(order)
    percluster.to_csv(f"{args.out}/{args.stem}_percluster.tsv", sep="\t")

    from scipy.stats import mannwhitneyu
    fc = str(args.focus_cluster)
    rows = []
    for col, label in [(basal_c, progs[0]), (emt_c, progs[1])]:
        a = df.loc[df["cluster"] == fc, col]
        b = df.loc[df["cluster"] != fc, col]
        u, p = mannwhitneyu(a, b, alternative="greater")
        rank = int((percluster[col].rank(ascending=False)[fc]))
        rows.append({"program": label, f"cl{fc}_mean": round(float(a.mean()), 4),
                     "rest_mean": round(float(b.mean()), 4),
                     "delta": round(float(a.mean() - b.mean()), 4),
                     "MWU_p_greater": f"{p:.3g}",
                     f"cl{fc}_rank_of_{len(order)}": rank})
    stats = pd.DataFrame(rows)
    stats.to_csv(f"{args.out}/{args.stem}_cl17_stats.tsv", sep="\t", index=False)
    print(f"[cl{fc} enrichment]\n" + stats.to_string(index=False))

    # ---- figure ----
    plt = setup_style(font_size=7)
    callouts = [g for g in args.callout_genes if g in df.columns]
    rowvars = [(basal_c, progs[0]), (emt_c, progs[1])] + [(g, g) for g in callouts]

    def draw_scatter(ax_sc):
        _draw_scatter(ax_sc, df, percluster, order, basal_c, emt_c, progs, colors, fc,
                      swap=args.scatter_swap_axes)

    def draw_violins(axes):
        _draw_violins(axes, df, order, rowvars, colors, fc)

    if args.panels == "scatter":
        # landscape default (was 5.8x5.4 square) -> fills the wide-short 2F cell; with --swap-axes
        # EMT (large range) on x uses the width. --figwidth-mm alone now yields a landscape panel.
        fig, ax_sc = plt.subplots(figsize=figsize(7.8, 4.6))
        draw_scatter(ax_sc)
    elif args.panels == "violins":
        n = len(rowvars)
        fig, axes = plt.subplots(n, 1, figsize=(max(6.0, 0.42 * len(order) + 2.0), 1.5 * n),
                                 sharex=True)
        draw_violins(np.atleast_1d(axes))
    elif args.panels == "scatter+tp63":
        tp = next((rv for rv in rowvars if rv[0] == "TP63"), None)
        if tp is None:
            print("[warn] TP63 not in object — rendering scatter only")
            fig, ax_sc = plt.subplots(figsize=(5.8, 5.4))
            draw_scatter(ax_sc)
        else:
            fig = plt.figure(figsize=(5.8, 6.6))
            gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.32)
            draw_scatter(fig.add_subplot(gs[0, 0]))
            _draw_violins([fig.add_subplot(gs[1, 0])], df, order, [tp], colors, fc)
    else:  # both
        n_right = len(rowvars)
        fig = plt.figure(figsize=(11, 6.5))
        gs = fig.add_gridspec(n_right, 2, width_ratios=[1.15, 1.0], wspace=0.28, hspace=0.55)
        draw_scatter(fig.add_subplot(gs[:, 0]))
        draw_violins([fig.add_subplot(gs[i, 1]) for i in range(n_right)])

    savefig(fig, args.out, args.stem, args.plot_format)
    plt.close(fig)
    print(f"[ok] wrote {args.stem}.{{{','.join(args.plot_format)}}} (panels={args.panels}) "
          f"+ 2 TSVs to {args.out}/")


def _draw_scatter(ax_sc, df, percluster, order, basal_c, emt_c, progs, colors, fc, swap=False):
    import numpy as np
    # swap -> EMT (larger range) on x, basal on y: fills a landscape cell without distortion
    xcol, ycol = (emt_c, basal_c) if swap else (basal_c, emt_c)
    xlab, ylab = (progs[1], progs[0]) if swap else (progs[0], progs[1])
    other = df["cluster"] != fc
    ax_sc.scatter(df.loc[other, xcol], df.loc[other, ycol], s=3, c="0.8",
                  alpha=0.35, edgecolor="none", zorder=1, rasterized=True)
    ax_sc.scatter(df.loc[~other, xcol], df.loc[~other, ycol], s=8,
                  c=colors.get(fc, "#d73027"), alpha=0.8, edgecolor="none",
                  zorder=3, label=f"cluster {fc}")
    # cluster centroids, labelled
    for c in order:
        m = percluster.loc[c]
        ax_sc.scatter(m[xcol], m[ycol], s=70, c=colors.get(c, "#444"),
                      edgecolor="black", lw=0.6, zorder=4)
        ax_sc.annotate(c, (m[xcol], m[ycol]), fontsize=6.5, ha="center", va="center",
                       zorder=5, color="white", fontweight="bold")
    ax_sc.axhline(0, color="0.6", lw=0.5, ls="--"); ax_sc.axvline(0, color="0.6", lw=0.5, ls="--")
    ax_sc.set_xlabel(f"{xlab} signature score")
    ax_sc.set_ylabel(f"{ylab} signature score")
    ax_sc.set_title(f"cluster {fc}: EMT-high / basal-intermediate \n— basal-like, not classical basal",
                    fontsize=6)
    # boxed legend (white bg + edge) so the single cluster-17 swatch reads as a legend,
    # not as another big dot on the scatter
    ax_sc.legend(loc="upper right", fontsize=7, frameon=True, facecolor="white",
                 edgecolor="0.7", framealpha=0.9, borderpad=0.4, handletextpad=0.4)
    for s in ("top", "right"):
        ax_sc.spines[s].set_visible(False)


def _draw_violins(axes, df, order, rowvars, colors, fc):
    import numpy as np
    x = np.arange(len(order))
    for i, (col, label) in enumerate(rowvars):
        ax = axes[i]
        data = [df.loc[df["cluster"] == c, col].values for c in order]
        parts = ax.violinplot(data, positions=x, widths=0.85, showextrema=False)
        for j, pc in enumerate(parts["bodies"]):
            c = order[j]
            pc.set_facecolor(colors.get(c, "#999"))
            pc.set_edgecolor("none")
            pc.set_alpha(0.95 if c == fc else 0.55)
        # highlight focus cluster position
        ax.axvline(order.index(fc), color="#d73027", lw=0.8, ls=":", zorder=0)
        ax.set_ylabel(label, fontsize=6.5)
        ax.set_xticks(x)
        if i == len(rowvars) - 1:
            ax.set_xticklabels(order, fontsize=6, rotation=90)
            if colors:
                for t, c in zip(ax.get_xticklabels(), order):
                    t.set_color(colors.get(c, "black"))
            ax.set_xlabel("cluster (adenocarcinoma → neuroendocrine)", fontsize=7)
        else:
            ax.set_xticklabels([])
        ax.tick_params(labelsize=6, length=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)


if __name__ == "__main__":
    main()
