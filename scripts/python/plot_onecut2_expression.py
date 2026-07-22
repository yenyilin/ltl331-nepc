#!/usr/bin/env python3
"""
plot_onecut2_expression.py — Fig 5 supplement: ONECUT2 DIRECT EXPRESSION.

Why this panel exists. ONECUT2 cannot be scored by either prior-based TF-activity method we
use: it was pruned by SCENIC's coexpression+motif construction (not in the
~200-regulon set) and it is ABSENT from the CollecTRI prior, so decoupleR returns
no activity estimate for it. We therefore answer the comment with direct expression:
ONECUT2 is shown per cluster alongside the canonical NEPC drivers that the
prior-based methods DO recover (ASCL1, NEUROD1, NKX2-1, POU3F2, FOXA2, INSM1), so
the reader sees ONECUT2 occupies the same neuroendocrine-cluster expression regime.

Dot SIZE  = fraction of cells in the cluster expressing the gene (>0).
Dot COLOR = mean expression. --scale none (DEFAULT) = raw mean lognorm (sequential
            Reds; the honest "is it expressed" view); --scale zscore = z across
            clusters (diverging RdBu_r), to match Fig 5A/5B colour convention.

ONECUT2 is highlighted (bold, leading row). Genes on y, clusters on x in the
canonical adeno->NE order; cluster tick labels coloured by the 18-shade palette.

CRITICAL: the base object's .X is SCALED — expression MUST come from .raw (lognorm).
Uses .raw by default. Reuses compute()/scale_cols() from plot_marker_dotplot.py.

Example:
  python plot_onecut2_expression.py --h5ad ltl331_base.h5ad \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --outdir figures --stem figS_onecut2_expression
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_marker_dotplot import compute, scale_cols  # noqa: E402  (sibling reuse, .raw-correct)

# ONECUT2 first (the R1#4 focus), then the canonical NEPC drivers the prior-based
# methods recover, as the comparison context.
DEFAULT_GENES = ['ONECUT2', 'ASCL1', 'NEUROD1', 'NKX2-1', 'POU3F2', 'FOXA2', 'INSM1']
FOCUS = 'ONECUT2'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--genes', nargs='+', default=DEFAULT_GENES,
                    help='genes to show (first one is highlighted as the focus)')
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--cluster-colors', default='data/cluster_colors_18.json',
                    help='cluster->hex json to colour the cluster tick labels')
    ap.add_argument('--use-raw', dest='use_raw', action='store_true', default=True)
    ap.add_argument('--no-use-raw', dest='use_raw', action='store_false')
    ap.add_argument('--scale', choices=['none', 'zscore', 'var'], default='none',
                    help='none = raw mean lognorm (default); zscore matches Fig 5A/5B')
    ap.add_argument('--vlim', type=float, default=2.0, help='z-score colour limit')
    ap.add_argument('--smin', type=float, default=4.0)
    ap.add_argument('--smax', type=float, default=200.0)
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='cluster ids to drop (default 18 = inferCNV reference)')
    ap.add_argument('--figwidth', type=float, default=None)
    ap.add_argument('--figheight', type=float, default=None)
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='figS_onecut2_expression')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm, Normalize
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            print(f'excluding clusters {args.exclude_clusters}: dropping {int((~keep).sum())} cells')
            adata = adata[keep.values].copy()

    # cluster order along the canonical axis
    present = sorted(adata.obs[args.cluster_key].astype(str).unique(),
                     key=lambda s: int(s) if s.isdigit() else 1e9)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        order = [c for c in co if c in present] + [c for c in present if c not in co]
    else:
        order = present

    mean_df, frac_df, missing = compute(adata, args.genes, args.cluster_key, order, args.use_raw)
    if missing:
        print(f'[warn] genes not found in object: {missing}')
    genes = [g for g in args.genes if g in mean_df.columns]
    if not genes:
        raise SystemExit('none of the requested genes were found in the object')
    if FOCUS not in genes and any(g.upper() == FOCUS for g in args.genes):
        print(f'[warn] {FOCUS} requested but NOT in the object var_names — '
              f'cannot show direct expression; check the gene symbol/alias')
    mean_df = mean_df[genes]; frac_df = frac_df[genes]

    # write the per-cluster stats for the response letter / supplement table
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    stats = pd.concat({'mean_lognorm': mean_df.T, 'frac_expressing': frac_df.T}, axis=1)
    stats.to_csv(outdir / f'{args.stem}_stats.tsv', sep='\t')
    if FOCUS in mean_df.columns:
        top = mean_df[FOCUS].astype(float).sort_values(ascending=False)
        print(f'[{FOCUS}] top clusters by mean lognorm: '
              + ', '.join(f'cl{c}={v:.2f}' for c, v in top.head(5).items()))

    color_df = scale_cols(mean_df, args.scale)
    if args.scale == 'zscore':
        cmap = plt.get_cmap('RdBu_r'); norm = TwoSlopeNorm(vmin=-args.vlim, vcenter=0, vmax=args.vlim)
        cbar_label = 'mean expression (z-score across clusters)'
    else:
        cmap = plt.get_cmap('Reds')
        norm = Normalize(float(np.nanmin(color_df.values)), float(np.nanmax(color_df.values)))
        cbar_label = ('relative expression (0–1 per gene)' if args.scale == 'var'
                      else 'mean expression (lognorm)')

    # genes on y (focus at top), clusters on x
    genes_y = list(genes)[::-1]  # so the first/focus gene plots at the TOP row
    ncl, ng = len(order), len(genes_y)
    fw = args.figwidth or max(5.0, 0.42 * ncl + 2.4)
    fh = args.figheight or max(2.4, 0.5 * ng + 1.6)
    fig, ax = plt.subplots(figsize=(fw, fh))

    palette = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if args.cluster_colors else {}
    for yi, gene in enumerate(genes_y):
        for xi, c in enumerate(order):
            f = frac_df.loc[c, gene]; v = color_df.loc[c, gene]
            if not np.isfinite(f) or not np.isfinite(v):
                continue
            ax.scatter(xi, yi, s=args.smin + (args.smax - args.smin) * float(f),
                       c=[cmap(norm(float(v)))], edgecolor='#444444', linewidth=0.3, zorder=2)

    ax.set_xticks(range(ncl))
    ax.set_xticklabels(order, fontsize=8)
    for tick, c in zip(ax.get_xticklabels(), order):
        tick.set_color(palette.get(c, '#000000'))
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes_y, fontsize=9, style='italic')
    for tick, g in zip(ax.get_yticklabels(), genes_y):
        if g.upper() == FOCUS:
            tick.set_fontweight('bold'); tick.set_style('italic')
    ax.set_xlim(-0.6, ncl - 0.4); ax.set_ylim(-0.6, ng - 0.4)
    ax.set_xlabel('cluster (adenocarcinoma → neuroendocrine axis)', fontsize=9)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title(f'Direct expression of canonical NEPC TFs — {FOCUS} highlighted\n'
                 f'({FOCUS} not scored by SCENIC/CollecTRI; shown by expression)',
                 fontsize=9)

    # colour bar + size legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(cbar_label, fontsize=7); cb.ax.tick_params(labelsize=6)
    for frac in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=args.smin + (args.smax - args.smin) * frac,
                   c='#888888', edgecolor='#444444', linewidth=0.3, label=f'{int(frac*100)}%')
    ax.legend(title='% expressing', loc='center left', bbox_to_anchor=(1.06, 0.5),
              fontsize=6, title_fontsize=7, frameon=False, labelspacing=1.1)

    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}} + {args.stem}_stats.tsv')


if __name__ == '__main__':
    main()
