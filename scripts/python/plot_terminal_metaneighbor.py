#!/usr/bin/env python3
"""
plot_terminal_metaneighbor.py — Fig 6B centerpiece: cross-cohort replicability of the
terminal/lineage states (MetaNeighbor-style neighbor-voting AUROC), EMBEDDING-FREE.

Two panels:
  left  — AUROC matrix (train group x test group), grouped by lineage then cohort,
          same-cohort cells blank (training must come from a different cohort).
  right — mean cross-cohort same-type AUROC per lineage, min-max whiskers, with the
          0.5 (chance) and 0.8 (replicable) reference lines.

Reads the TSVs written by terminal_metaneighbor.py. No object needed.

Example:
  python scripts/plot_terminal_metaneighbor.py \
      --indir figures/data/terminal_metaneighbor --outdir figures/fig6
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# display order: lineage (adeno -> NE terminals -> basal/other), cohort within
CT_ORDER = ['PRAD', 'NE_ASCL1pos', 'NE_ASCL1neg', 'Basal', 'Other']
COHORT_ORDER = ['LTL331', 'Gao', 'Li']
PRETTY_CT = {'PRAD': 'Adeno', 'NE_ASCL1pos': 'NE ASCL1+', 'NE_ASCL1neg': 'NE ASCL1−',
             'Basal': 'Basal', 'Other': 'Other'}


def order_groups(groups):
    def key(g):
        coh, ct = g.split('|', 1)
        return (CT_ORDER.index(ct) if ct in CT_ORDER else 99,
                COHORT_ORDER.index(coh) if coh in COHORT_ORDER else 99)
    return sorted(groups, key=key)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--indir', default='figures/data/terminal_metaneighbor')
    p.add_argument('--outdir', default='figures/fig6')
    p.add_argument('--stem', default='fig6B_metaneighbor')
    p.add_argument('--figwidth-mm', type=float, default=None,
                   help='author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa')
    p.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm')
    args = p.parse_args()

    indir, outdir = Path(args.indir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

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

    mat = pd.read_csv(indir / 'metaneighbor_auroc_matrix.tsv', sep='\t', index_col=0)
    summ = pd.read_csv(indir / 'metaneighbor_crossdataset_summary.tsv', sep='\t')

    order = order_groups(list(mat.index))
    mat = mat.reindex(index=order, columns=order)

    from _utils import setup_style, savefig
    plt = setup_style()
    import matplotlib.pyplot as mpl
    fig, (axH, axB) = plt.subplots(1, 2, figsize=figsize(9.4, 4.4),
                                   gridspec_kw={'width_ratios': [1.55, 1]})

    # ---- left: AUROC heatmap ------------------------------------------------
    M = mat.values.astype(float)
    im = axH.imshow(np.ma.masked_invalid(M), cmap='RdBu_r', vmin=0.0, vmax=1.0,
                    aspect='equal')
    labels = [f'{PRETTY_CT.get(g.split("|")[1], g.split("|")[1])}\n{g.split("|")[0]}'
              for g in order]
    axH.set_xticks(range(len(order))); axH.set_xticklabels(labels, fontsize=3.5, rotation=90)
    axH.set_yticks(range(len(order))); axH.set_yticklabels(labels, fontsize=3.5)
    axH.set_xlabel('test group (cohort | state)', fontsize=6.5)
    axH.set_ylabel('train group', fontsize=6.5)
    axH.set_title('neighbor-voting AUROC', fontsize=10)#16)
    for i in range(len(order)):
        for j in range(len(order)):
            v = M[i, j]
            if not np.isnan(v):
                axH.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=3.5,
                         color='white' if (v > 0.78 or v < 0.30) else 'black')
    # lineage block separators
    seps = np.cumsum([sum(1 for g in order if g.split('|')[1] == ct) for ct in CT_ORDER])[:-1]
    for s in seps:
        axH.axhline(s - 0.5, color='0.3', lw=0.5); axH.axvline(s - 0.5, color='0.3', lw=0.5)
    cb = fig.colorbar(im, ax=axH, fraction=0.046, pad=0.02)
    cb.set_label('AUROC', fontsize=5); cb.ax.tick_params(labelsize=4.5)

    # ---- right: summary bar -------------------------------------------------
    summ['ct'] = pd.Categorical(summ['cell_type'], CT_ORDER, ordered=True)
    summ = summ.sort_values('ct')
    y = np.arange(len(summ))[::-1]
    colors = ['#2c7fb8' if v >= 0.8 else ('#9ecae1' if v >= 0.65 else '#cccccc')
              for v in summ['mean_cross_dataset_AUROC']]
    axB.barh(y, summ['mean_cross_dataset_AUROC'], color=colors, height=0.62,
             xerr=[summ['mean_cross_dataset_AUROC'] - summ['min'],
                   summ['max'] - summ['mean_cross_dataset_AUROC']],
             error_kw=dict(lw=0.7, capsize=2, ecolor='0.35'))
    axB.axvline(0.5, color='0.5', lw=0.7, ls=':')
    axB.axvline(0.8, color='#d95f02', lw=0.8, ls='--')
    axB.text(0.5, -0.7, 'chance', fontsize=4.5, color='0.5', ha='center')
    axB.text(0.8, -0.7, 'replicable', fontsize=4.5, color='#d95f02', ha='center')
    for yi, v, mx in zip(y, summ['mean_cross_dataset_AUROC'], summ['max']):
        axB.text(mx + 0.03, yi, f'{v:.2f}', va='center', fontsize=4.5)
    axB.set_yticks(y)
    axB.set_yticklabels([PRETTY_CT.get(c, c) for c in summ['cell_type']], fontsize=6.5)#8.5)
    axB.set_xlim(0, 1.05); axB.set_xlabel('mean cross-cohort AUROC', fontsize=6.5)#9)
    axB.set_title('cross-cohort cell-state replicability', fontsize=10)#16)
    axB.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    savefig(fig, outdir, args.stem, ('pdf', 'png'))
    plt.close(fig)
    print(f'[6B] ASCL1+ NE={summ.set_index("cell_type").loc["NE_ASCL1pos","mean_cross_dataset_AUROC"]:.3f} '
          f'PRAD={summ.set_index("cell_type").loc["PRAD","mean_cross_dataset_AUROC"]:.3f} '
          f'ASCL1− NE={summ.set_index("cell_type").loc["NE_ASCL1neg","mean_cross_dataset_AUROC"]:.3f}')


if __name__ == '__main__':
    main()
