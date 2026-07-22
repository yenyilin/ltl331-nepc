#!/usr/bin/env python3
"""
plot_timepoint_panel.py — the Fig 2 temporal-axis panel(s): how the cancer-cell
population is distributed across the 8 NEtD timepoints. DESCRIPTIVE only (no arrows,
no pseudotime) — directionality belongs to the trajectory figure.

  --which umap         F1: faceted UMAP, one mini-plot per timepoint, that timepoint's
                       cells highlighted over a gray background of all cells. The literal
                       "timepoint coloring", done without 8-colour overplotting; shows the
                       cloud traversing transcriptional space (AR/adeno -> NE) over time.
  --which composition  F2: per-timepoint cluster-composition stacked bar (each bar = a
                       timepoint, sums to 1), using the SAME cluster colours as Fig 2A.
                       Quantitative; the mid-course bulge of intermediate clusters (4/6/17)
                       is the visual evidence of transition through an intermediate.
  --which both         both (default).

Timepoints are ordered CHRONOLOGICALLY by week_num (so wk2 precedes wk12, etc.).
Cluster stacking order follows data/cluster_order.json; cluster colours follow
data/cluster_phenotype_colors.json (falls back to tab20). Reads obs only + obsm UMAP —
no expression matrix needed, so the scaled-.X issue does not apply here.

Example:
  python plot_timepoint_panel.py --h5ad ltl331_base.h5ad \\
      --cluster-key seurat_clusters --timepoint-key timepoint --week-key week_num \\
      --cluster-order data/cluster_order.json --colors data/cluster_phenotype_colors.json \\
      --which both --outdir figures --stem fig2_timepoint
"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd

# coarse phenotype bands for --by-phenotype (ordered adeno -> NE), with pre-EMT split out
# so the transition reads fully. Colours extend data/cluster_phenotype_colors.json
# (pre-EMT gets its own amber since that file lumps it into AR-low orange).
PHENO_ORDER = ['AR-high PRAD', 'AR-low PRAD', 'Pre-EMT', 'Intermediate (EMT)', 'NEPC']
PHENO_COLOR = {'AR-high PRAD': '#7570b3', 'AR-low PRAD': '#d95f02', 'Pre-EMT': '#e6ab02',
               'Intermediate (EMT)': '#e7298a', 'NEPC': '#1b9e77'}


def phenotype_of(label):
    l = label.lower()
    if 'ar-high' in l:       return 'AR-high PRAD'
    if 'pre-emt' in l:       return 'Pre-EMT'
    if 'ar-low' in l:        return 'AR-low PRAD'
    if 'intermediate' in l:  return 'Intermediate (EMT)'
    if 'nepc' in l or 'ne ' in l: return 'NEPC'
    return 'other'


def order_timepoints(obs, tp_key, week_key):
    tp = obs[tp_key].astype(str)
    if week_key in obs:
        wk = pd.to_numeric(obs[week_key], errors='coerce').groupby(tp).median()
        return list(wk.sort_values().index), wk
    # natural sort: preCx first, then by embedded digits
    def k(s):
        if 'pre' in s.lower(): return -1
        d = ''.join(c for c in s if c.isdigit())
        return int(d) if d else 999
    return sorted(tp.unique(), key=k), None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--timepoint-key', default='timepoint')
    ap.add_argument('--week-key', default='week_num')
    ap.add_argument('--umap-key', default='X_umap')
    ap.add_argument('--cluster-order', default=None)
    ap.add_argument('--colors', default=None, help='cluster_phenotype_colors.json')
    ap.add_argument('--by-phenotype', action='store_true',
                    help='F2: collapse clusters to coarse phenotype bands (clean adeno->NE shift)')
    ap.add_argument('--labels', default='data/cluster_phenotype_labels.json',
                    help='cluster->phenotype labels json (for --by-phenotype)')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='cluster ids to drop (default 18 = inferCNV reference)')
    ap.add_argument('--which', choices=['umap', 'composition', 'both'], default='both')
    ap.add_argument('--cmap', default='viridis',
                    help='sequential matplotlib colormap for the timepoint UMAP gradient '
                         '(e.g. viridis, plasma, cividis). Ordered data -> keep it sequential.')
    ap.add_argument('--cmap-hi', type=float, default=1.0,
                    help='upper fraction of the colormap to use (0<hi<=1); e.g. 0.72 ends the '
                         'ramp on green instead of yellow so the last timepoint echoes the NEPC '
                         'cluster colour. Default 1.0 = full range.')
    ap.add_argument('--ncols', type=int, default=4, help='facets per row for the UMAP panel')
    ap.add_argument('--bg-cap', type=int, default=30000,
                    help='max gray background cells per facet (render speed)')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig2_timepoint')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author panel at this FINAL width in mm (assemble_figure --qa onpage_mm '
                         'target) so it scales ~1:1; the saved PDF is wider by any outside legend — '
                         'tune with --qa. This is the TOTAL panel width (facets divide it).')
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
    tp_order, wk = order_timepoints(obs, args.timepoint_key, args.week_key)
    print(f'timepoint order (chronological): {tp_order}')

    # cluster stacking order + colours
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

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # ---- F1: faceted UMAP ----
    if args.which in ('umap', 'both'):
        if args.umap_key not in adata.obsm:
            raise SystemExit(f'{args.umap_key} not in adata.obsm ({list(adata.obsm)})')
        um = np.asarray(adata.obsm[args.umap_key])[:, :2]
        n = len(tp_order)
        ncols = min(args.ncols, n)
        nrows = int(np.ceil(n / ncols))
        hi = min(max(args.cmap_hi, 1e-3), 1.0)
        if hi != args.cmap_hi:
            print(f'[warn] --cmap-hi {args.cmap_hi} clamped to {hi} (expected 0<hi<=1)')
        vir = plt.get_cmap(args.cmap)
        tp_color = {t: vir((i / max(n - 1, 1)) * hi) for i, t in enumerate(tp_order)}
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=figsize(2.6 * ncols, 2.6 * nrows),
                                 squeeze=False)
        bg = um
        if len(bg) > args.bg_cap:
            bg_idx = rng.choice(len(bg), args.bg_cap, replace=False)
        else:
            bg_idx = np.arange(len(bg))
        for i, t in enumerate(tp_order):
            ax = axes[i // ncols][i % ncols]
            ax.scatter(um[bg_idx, 0], um[bg_idx, 1], s=1.5, c='#dddddd',
                       lw=0, rasterized=True, zorder=1)
            m = tp_all == t
            ax.scatter(um[m, 0], um[m, 1], s=2.5, color=tp_color[t],
                       lw=0, rasterized=True, zorder=2)
            ax.set_title(f'{t}  (n={int(m.sum())})', fontsize=5.5, pad=2)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
        for j in range(n, nrows * ncols):
            axes[j // ncols][j % ncols].axis('off')
        fig.suptitle('Cell distribution across NEtD timepoints (UMAP)', fontsize=8, y=0.99)
        for fmt in args.plot_format:
            fig.savefig(outdir / f'{args.stem}_umap.{fmt}', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'[ok] wrote {args.stem}_umap.{{{",".join(args.plot_format)}}}')

    # ---- F2: per-timepoint composition (each bar sums to 1) ----
    if args.which in ('composition', 'both'):
        if args.by_phenotype:
            labels = {str(k): v for k, v in json.load(open(args.labels)).items()}
            stack_vals = np.array([phenotype_of(labels.get(c, '')) for c in clusters])
            unknown = sorted(set(clusters[stack_vals == 'other']))
            if unknown:
                print(f'[warn] clusters with no phenotype label (shown as "other"): {unknown}')
            seen = set(stack_vals)
            stack_order = [p for p in PHENO_ORDER if p in seen] + (['other'] if 'other' in seen else [])
            stack_colors = dict(PHENO_COLOR, other='#cccccc')
            legend_title, title = 'phenotype', 'Phenotype composition per timepoint'
        else:
            stack_vals = clusters
            stack_order = cl_order
            stack_colors = cmap_cl
            legend_title, title = 'cluster', 'Cluster composition per timepoint'
        ct = pd.crosstab(pd.Series(tp_all, name='tp'),
                         pd.Series(stack_vals, name='grp'), normalize='index')
        ct = ct.reindex(index=tp_order, columns=stack_order, fill_value=0)
        fig, ax = plt.subplots(figsize=figsize(max(5, 0.7 * len(tp_order) + 2), 4.2))
        bottom = np.zeros(len(tp_order))
        for k in stack_order:
            ax.bar(range(len(tp_order)), ct[k].values, bottom=bottom, width=0.82,
                   color=stack_colors.get(k, '#cccccc'), edgecolor='white', lw=0.2, label=k)
            bottom += ct[k].values
        ax.set_xticks(range(len(tp_order)))
        ax.set_xticklabels(tp_order, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('fraction of cells', fontsize=16)
        ax.set_ylim(0, 1); ax.set_xlim(-0.6, len(tp_order) - 0.4)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.set_title(title, fontsize=24)
        handles = [Patch(facecolor=stack_colors.get(k, '#cccccc'), label=k) for k in stack_order]
        ax.legend(handles=handles, title=legend_title, loc='center left',
                  bbox_to_anchor=(1.005, 0.5), fontsize=16,
                  title_fontsize=16, ncol=1, frameon=False)
        for fmt in args.plot_format:
            fig.savefig(outdir / f'{args.stem}_composition.{fmt}', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'[ok] wrote {args.stem}_composition.{{{",".join(args.plot_format)}}}')


if __name__ == '__main__':
    main()
