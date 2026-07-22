#!/usr/bin/env python3
"""
fig3d_nes_heatmap.py — publication Fig 3D: per-cluster Hallmark GSEA NES heatmap,
columns ordered along the NEtD trajectory, rows ordered by PRAD->NEPC trend, with a
phenotype colour strip, coarse group headers, and FDR significance dots.

Rank-based NES (from per_cluster_gsea.py) — not subject to the highly-expressed-gene
bias of per-cell module scores. Drop-in for the manuscript Fig 3D.

Inputs (all from per_cluster_gsea.py + repo data/):
  per_cluster_NES.tsv, per_cluster_FDR.tsv, cluster_phenotype_labels.json,
  cluster_phenotype_colors.json, per_cell_trajectory.tsv (for velocity_pt order).

Example:
  python fig3d_nes_heatmap.py \
    --nes data/per_cluster_gsea/per_cluster_NES.tsv \
    --fdr data/per_cluster_gsea/per_cluster_FDR.tsv \
    --labels data/cluster_phenotype_labels.json \
    --colors data/cluster_phenotype_colors.json \
    --order-by data/scanpy/trajectory_summary_3term/per_cell_trajectory.tsv \
    --outdir data/per_cluster_gsea --plot-format pdf png
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd


def group_of(label):
    l = label.lower()
    if 'ar-high prad' in l: return 'AR-high PRAD'
    if 'pre-emt' in l:      return 'AR-low PRAD (pre-EMT)'
    if 'ar-low prad' in l:  return 'AR-low PRAD'
    if 'intermediate' in l: return 'Intermediate'
    if 'nepc' in l:         return 'NEPC'
    return 'other'


FINE_ORDER = ['AR-high PRAD', 'AR-low PRAD', 'AR-low PRAD (pre-EMT)', 'Intermediate', 'NEPC']
COARSE = {'AR-high PRAD': 'PRAD', 'AR-low PRAD': 'PRAD', 'AR-low PRAD (pre-EMT)': 'PRAD',
          'Intermediate': 'Intermediate', 'NEPC': 'NEPC'}
COARSE_ORDER = ['PRAD', 'Intermediate', 'NEPC']


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nes', required=True)
    ap.add_argument('--fdr', default=None)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--colors', default=None)
    ap.add_argument('--cluster-order', default=None,
                    help='cluster_order.json with an "order" list (phenotypic axis). '
                         'Manuscript default — overrides --order-by pseudotime sorting.')
    ap.add_argument('--order-by', default=None,
                    help='per_cell_trajectory.tsv (pseudotime order; for Fig 4-style only, '
                         'not the manuscript-wide axis)')
    ap.add_argument('--order-col', default='velocity_pt')
    ap.add_argument('--order-cluster-col', default='cluster')
    ap.add_argument('--fdr-sig', type=float, default=0.05, help='dot if FDR < this')
    ap.add_argument('--vmax', type=float, default=None)
    ap.add_argument('--highlight', nargs='+',
                    default=['KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN'],
                    help='pathway substrings to bold on the y-axis')
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--name', default='fig3D_gsea_NES')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    nes = pd.read_csv(args.nes, sep='\t', index_col=0); nes.columns = nes.columns.map(str)
    fdr = pd.read_csv(args.fdr, sep='\t', index_col=0) if args.fdr else None
    if fdr is not None: fdr.columns = fdr.columns.map(str)
    labels = {str(k): v for k, v in json.load(open(args.labels)).items()}
    colors = {str(k): v for k, v in json.load(open(args.colors)).items()} if args.colors else {}
    outdir = Path(args.outdir) if args.outdir else Path(args.nes).parent
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- column order ----
    cl_fine = {c: group_of(labels.get(c, '')) for c in nes.columns}
    xlab = 'cluster'
    if args.cluster_order:                       # manuscript default: canonical phenotypic axis
        canon = json.load(open(args.cluster_order))
        order = [str(c) for c in canon.get('order', [])]
        cols = [c for c in order if c in nes.columns]
        extra = [c for c in nes.columns if c not in cols and cl_fine[c] in FINE_ORDER]
        if extra:
            print(f'[warn] clusters absent from cluster_order.json, appended: {extra}')
            cols += extra
        xlab = canon.get('axis_label', 'cluster (adenocarcinoma → neuroendocrine phenotypic axis)')
    else:                                        # legacy: phenotype group, pseudotime within
        cols = [c for c in nes.columns if cl_fine[c] in FINE_ORDER]
        if args.order_by:
            od = pd.read_csv(args.order_by, sep='\t')
            med = od.groupby(args.order_cluster_col)[args.order_col].median()
            med.index = med.index.map(str)
            okey = {c: float(med.get(c, np.inf)) for c in cols}
            xlab = 'cluster (ordered by NEtD stage; pseudotime within stage)'
        else:
            okey = {c: (int(c) if c.isdigit() else 1e9) for c in cols}
        cols.sort(key=lambda c: (FINE_ORDER.index(cl_fine[c]), okey[c]))

    # ---- row order: PRAD->NEPC trend (NEPC mean - PRAD mean), descending ----
    coarse = {c: COARSE[cl_fine[c]] for c in cols}
    gm = nes[cols].T.groupby(pd.Series(coarse)).mean().T
    trend = (gm.get('NEPC', 0) - gm.get('PRAD', 0)).sort_values(ascending=False)
    rows = list(trend.index)

    M = nes.loc[rows, cols]
    F = fdr.loc[rows, cols] if fdr is not None else None
    vmax = args.vmax or float(np.ceil(np.nanmax(np.abs(M.values)) * 2) / 2)
    vmax = min(vmax, 3.5)

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42   # editable text in Illustrator
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    nrow, ncol = M.shape
    W = 0.42 * ncol + 5.5
    H = 0.22 * nrow + 2.8
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, nrow * 0.22],
                          width_ratios=[0.42 * ncol, 0.5], hspace=0.04, wspace=0.04)
    ax_strip = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])   # no sharex: hiding strip ticks must not clear heatmap labels
    ax_cb = fig.add_subplot(gs[1, 1])

    # heatmap
    im = ax.imshow(M.values, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(ncol)); ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(nrow))
    ylabels = [r.replace('HALLMARK_', '').replace('_', ' ').title() for r in rows]
    ax.set_yticklabels(ylabels, fontsize=7)
    for t, r in zip(ax.get_yticklabels(), rows):
        if any(h in r for h in args.highlight):
            t.set_fontweight('bold')
    ax.set_xlabel(xlab, fontsize=9)

    # FDR significance dots
    if F is not None:
        ys, xs = np.where(F.values < args.fdr_sig)
        ax.scatter(xs, ys, s=6, c='k', marker='o', linewidths=0)

    # coarse-group separators
    seps = [i for i in range(1, ncol) if coarse[cols[i]] != coarse[cols[i - 1]]]
    for s in seps:
        ax.axvline(s - 0.5, color='k', lw=1.3)

    # phenotype colour strip + coarse headers
    strip = np.ones((1, ncol, 3))
    for j, c in enumerate(cols):
        hexc = colors.get(c, '#cccccc').lstrip('#')
        strip[0, j] = [int(hexc[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    # rasterize the RGBA strip so the PDF embeds a raster of it (matches the PNG);
    # interpolation='nearest' keeps cell edges crisp and dodges the PDF-blank-band bug
    ax_strip.imshow(strip, aspect='auto', interpolation='nearest', rasterized=True)
    ax_strip.set_yticks([]); ax_strip.set_xticks([])
    ax_strip.set_xlim(-0.5, ncol - 0.5)   # align strip to heatmap columns (no sharex)
    ax_strip.set_ylim(0.5, -1.6)
    for sp in ax_strip.spines.values():
        sp.set_visible(False)
    # coarse group header text above the strip
    for g in COARSE_ORDER:
        idx = [j for j, c in enumerate(cols) if coarse[c] == g]
        if not idx: continue
        ax_strip.text(np.mean(idx), -1.0, g, ha='center', va='center',
                      fontsize=11, fontweight='bold')

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label('GSEA NES', fontsize=9)

    fig.suptitle('Hallmark pathway enrichment across the NEtD trajectory '
                 f'(rank-based GSEA per cluster; dot = FDR < {args.fdr_sig})',
                 fontsize=11, y=0.998)

    # phenotype legend at the bottom (below the x-axis label)
    seen, handles = set(), []
    for c in cols:
        f = cl_fine[c]
        if f in seen: continue
        seen.add(f)
        handles.append(Patch(facecolor=colors.get(c, '#cccccc'), label=f))
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.02),
               ncol=len(handles), fontsize=8, frameon=False)

    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.name}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.name}.{{{",".join(args.plot_format)}}} to {outdir}/'
          f'  ({nrow} pathways x {ncol} clusters, vmax={vmax})')
    print('col order:', ' '.join(cols))


if __name__ == '__main__':
    main()
