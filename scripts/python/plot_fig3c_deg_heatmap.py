#!/usr/bin/env python3
"""
plot_fig3c_deg_heatmap.py - Top-N DEG heatmap per cluster (original Fig 3C).
Multi-pairwise (1-vs-rest) rank_genes_groups with Wilcoxon, matching Seurat's
FindAllMarkers default. Picks top-N markers per cluster, deduplicates while
preserving cluster order, plots a z-scored expression heatmap. Optionally
bold/color-highlights a curated gene list via --highlight.

Outputs (in --outdir):
  - fig3c_deg_heatmap.{pdf,png}
  - fig3c_top{N}_per_cluster.tsv
  - fig3c_rank_genes_groups_full.tsv
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import scanpy as sc

from _utils import setup_style, savefig, numkey


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, help="Input AnnData file.")
    p.add_argument("--outdir", default="scripts/fig3c")
    p.add_argument("--cluster-col", default="seurat_clusters",
                   help="obs column with cluster labels.")
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("--method", default="wilcoxon",
                   choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"])
    p.add_argument("--use-raw", action="store_true",
                   help="Use adata.raw for rank_genes_groups (matches Seurat FindAllMarkers).")
    p.add_argument("--highlight", nargs="*", default=[],
                   help="Genes to bold/color on the y-axis (case-sensitive).")
    p.add_argument("--cluster-order", nargs="*", default=None,
                   help="Explicit cluster order; default = numeric ascending.")
    p.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    return p.parse_args()


def main():
    args = parse_args()
    plt = setup_style(font_size=8)
    os.makedirs(args.outdir, exist_ok=True)

    print(f"reading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    if args.cluster_col not in adata.obs.columns:
        raise SystemExit(
            f"--cluster-col '{args.cluster_col}' not in obs. Available: "
            + ", ".join(sorted(adata.obs.columns))
        )

    if args.cluster_order is not None:
        order = list(args.cluster_order)
    else:
        order = sorted(adata.obs[args.cluster_col].astype(str).unique(), key=numkey)
    adata.obs[args.cluster_col] = pd.Categorical(
        adata.obs[args.cluster_col].astype(str),
        categories=order, ordered=True,
    )

    print(f"rank_genes_groups (method={args.method}, use_raw={args.use_raw})")
    sc.tl.rank_genes_groups(
        adata, groupby=args.cluster_col, method=args.method,
        n_genes=200, use_raw=args.use_raw, pts=True,
    )

    full = sc.get.rank_genes_groups_df(adata, group=None)
    full_out = os.path.join(args.outdir, "fig3c_rank_genes_groups_full.tsv")
    full.to_csv(full_out, sep="\t", index=False)
    print(f"  wrote {full_out}")

    top = (full.groupby("group", sort=False)
                .head(args.top_n)
                .reset_index(drop=True))
    top_out = os.path.join(args.outdir, f"fig3c_top{args.top_n}_per_cluster.tsv")
    top.to_csv(top_out, sep="\t", index=False)
    print(f"  wrote {top_out}")

    gene_order, seen = [], set()
    for c in order:
        for g in top.loc[top["group"].astype(str) == str(c), "names"].tolist():
            if g not in seen and g in adata.var_names:
                gene_order.append(g)
                seen.add(g)
    print(f"plotting heatmap: {len(gene_order)} unique genes × {len(order)} clusters")

    sc.pl.heatmap(
        adata, var_names=gene_order, groupby=args.cluster_col,
        use_raw=args.use_raw, standard_scale="var",
        swap_axes=True, show_gene_labels=True, dendrogram=False,
        figsize=(max(6.0, len(order) * 0.45),
                 max(8.0, len(gene_order) * 0.16)),
        show=False,
    )
    fig = plt.gcf()

    if args.highlight:
        hi = set(args.highlight)
        for ax in fig.axes:
            for lbl in ax.get_yticklabels():
                if lbl.get_text() in hi:
                    lbl.set_fontweight("bold")
                    lbl.set_color("#B2182B")

    savefig(fig, args.outdir, "fig3c_deg_heatmap", formats=args.plot_format)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
