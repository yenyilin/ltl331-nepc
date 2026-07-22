#!/usr/bin/env python3
"""
plot_cluster_proportions.py — Fig 3A: cluster proportions across the 8 NEtD
timepoints. The "cluster dynamics" panel — how each transcriptional state rises
and falls over the adenocarcinoma -> NE course (the mid-course bulge of the
intermediate clusters 4/6/17 is the visual point of the figure).

Uses the SAME cluster colours (data/cluster_colors_18.json) and SAME stacking
order (data/cluster_order.json) as Fig 2A / the Fig 2-3 composition panel, so a
cluster's colour is identical across every figure.

  --mode stacked   each bar = one timepoint, clusters stacked (default). Pair with
                   --counts for absolute cell numbers instead of fractions.
  --mode area      stacked area across timepoints — smoother view of the same data,
                   emphasises the rise/fall dynamics of each cluster over time.

Timepoints are ordered CHRONOLOGICALLY by week_num. Reads obs only (no expression
matrix), so the scaled-.X issue does not apply.

Example:
  python plot_cluster_proportions.py --h5ad ltl331_base.h5ad \\
      --cluster-order data/cluster_order.json --colors data/cluster_colors_18.json \\
      --outdir figures --stem fig3A_cluster_proportions
"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd


def order_timepoints(obs, tp_key, week_key):
    """Chronological timepoint order via median week_num; natural-sort fallback."""
    tp = obs[tp_key].astype(str)
    if week_key in obs:
        wk = pd.to_numeric(obs[week_key], errors='coerce').groupby(tp).median()
        return list(wk.sort_values().index)
    def k(s):
        if 'pre' in s.lower():
            return -1
        d = ''.join(c for c in s if c.isdigit())
        return int(d) if d else 999
    return sorted(tp.unique(), key=k)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--timepoint-key', default='timepoint')
    ap.add_argument('--week-key', default='week_num')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--colors', default='data/cluster_colors_18.json',
                    help='flat {cluster: hex} palette (18-shade phenotype palette)')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='cluster ids to drop (default 18 = inferCNV reference)')
    ap.add_argument('--mode', choices=['stacked', 'area'], default='stacked')
    ap.add_argument('--counts', action='store_true',
                    help='plot absolute cell counts instead of per-timepoint fractions')
    ap.add_argument('--label-clusters', action='store_true',
                    help='write the cluster id inside each band at its thickest timepoint '
                         '(least distorted point); bands thinner than --label-min are skipped '
                         'and left to the legend. Works in both --mode area and stacked.')
    ap.add_argument('--label-min', type=float, default=0.04,
                    help='min band thickness to label, as a fraction of that timepoint\'s '
                         'cells (default 0.04 = 4%%). Raise to declutter, lower to label more.')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig3A_cluster_proportions')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa')
    ap.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm')
    args = ap.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh:
            return ((fw / 25.4) * w_mult, (fh / 25.4) * h_mult)
        if fw:   # width only -> preserve authored aspect ratio (avoids skew + empty margins)
            w = (fw / 25.4) * w_mult
            return (w, default_h_in * w / default_w_in)
        if fh:   # height only -> preserve authored aspect ratio
            h = (fh / 25.4) * h_mult
            return (default_w_in * h / default_h_in, h)
        return (default_w_in, default_h_in)

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            print(f'excluding clusters {args.exclude_clusters}: dropping {int((~keep).sum())} cells')
            adata = adata[keep.values].copy()
    obs = adata.obs
    clusters = obs[args.cluster_key].astype(str).values
    tp_all = obs[args.timepoint_key].astype(str).values
    tp_order = order_timepoints(obs, args.timepoint_key, args.week_key)
    print(f'timepoint order (chronological): {tp_order}')

    # cluster stacking order + colours (same source as every other figure)
    present = sorted(set(clusters), key=lambda s: int(s) if s.isdigit() else 1e9)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        cl_order = [c for c in co if c in present] + [c for c in present if c not in co]
    else:
        cl_order = present
    if args.colors:
        cmap_cl = {str(k): v for k, v in json.load(open(args.colors)).items()}
    else:
        tab = plt.get_cmap('tab20')
        cmap_cl = {c: tab(i % 20) for i, c in enumerate(cl_order)}
    missing = [c for c in cl_order if c not in cmap_cl]
    if missing:
        print(f'[warn] no colour for clusters {missing} -> grey')

    # timepoint x cluster table (counts, then optionally normalise per timepoint)
    ct = pd.crosstab(pd.Series(tp_all, name='tp'), pd.Series(clusters, name='cl'))
    ct = ct.reindex(index=tp_order, columns=cl_order, fill_value=0)
    if not args.counts:
        ct = ct.div(ct.sum(axis=1), axis=0).fillna(0)
    ylabel = 'number of cells' if args.counts else 'fraction of cells'

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    x = np.arange(len(tp_order))

    fig, ax = plt.subplots(figsize=figsize(max(5, 0.7 * len(tp_order) + 2), 4.2))
    # M[j] = cluster cl_order[j] across timepoints; same stacking order for both modes
    M = np.vstack([ct[c].values for c in cl_order])
    if args.mode == 'stacked':
        bottom = np.zeros(len(tp_order))
        for j, c in enumerate(cl_order):
            ax.bar(x, M[j], bottom=bottom, width=0.82,
                   color=cmap_cl.get(c, '#cccccc'), edgecolor='white', lw=0.2, label=c)
            bottom += M[j]
    else:  # area
        ax.stackplot(x, M, colors=[cmap_cl.get(c, '#cccccc') for c in cl_order],
                     labels=cl_order, edgecolor='white', linewidth=0.2)

    if args.label_clusters:
        import matplotlib.patheffects as pe
        cum = np.cumsum(M, axis=0)          # top edge of each band per timepoint
        coltot = M.sum(axis=0)              # per-timepoint total (=1 in fraction mode)
        for j, c in enumerate(cl_order):
            ti = int(np.argmax(M[j]))       # the band's thickest (least-distorted) timepoint
            thick = M[j, ti]
            if coltot[ti] <= 0 or thick / coltot[ti] < args.label_min:
                continue                    # too thin here -> leave it to the legend
            ax.text(x[ti], cum[j, ti] - thick / 2, c, ha='center', va='center',
                    fontsize=7, color='white', fontweight='bold', zorder=6,
                    path_effects=[pe.withStroke(linewidth=1.4, foreground='black')])

    ax.set_xticks(x)
    ax.set_xticklabels(tp_order, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xlim(-0.6, len(tp_order) - 0.4)
    if not args.counts:
        ax.set_ylim(0, 1)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title('Cluster proportions across\nNEtD timepoints', fontsize=11)
    handles = [Patch(facecolor=cmap_cl.get(c, '#cccccc'), label=c) for c in cl_order]
    ax.legend(handles=handles, title='cluster', loc='center left',
              bbox_to_anchor=(1.005, 0.5), fontsize=6.5, title_fontsize=8,
              ncol=1, frameon=False)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}_{args.mode}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}_{args.mode}.{{{",".join(args.plot_format)}}}')


if __name__ == '__main__':
    main()
