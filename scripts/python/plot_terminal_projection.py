#!/usr/bin/env python3
"""
plot_terminal_projection.py — Fig 6C: reference-projection label transfer. Clinical cells
projected onto the fixed LTL331 manifold inherit the LTL331 terminal labels.

Two panels:
  left  — confusion heatmap: query ORIGINAL coarse lineage (rows) x TRANSFERRED LTL331
          label (cols), row-normalized. Strong diagonal on PRAD and NE = label inheritance.
  right — per-cohort landing of query NE cells across transferred labels (stacked fraction).

Reads the TSVs written by terminal_reference_project.py. No h5ad object needed.

Example:
  python scripts/plot_terminal_projection.py \
      --indir figures/data/terminal_projection --outdir figures/fig6
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

COL_ORDER = ['PRAD', 'NE_ASCL1pos', 'NE_ASCL1neg', 'Basal', 'Other']
ROW_ORDER = ['PRAD', 'NE', 'Basal', 'Other']
PRETTY = {'PRAD': 'Adeno', 'NE': 'NE', 'NE_ASCL1pos': 'NE ASCL1+',
          'NE_ASCL1neg': 'NE ASCL1−', 'Basal': 'Basal', 'Other': 'Other'}
TCOLORS = {'PRAD': '#4575b4', 'NE_ASCL1pos': '#d73027', 'NE_ASCL1neg': '#fc8d59',
           'Basal': '#91bfdb', 'Other': '#cccccc'}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--indir', default='figures/data/terminal_projection')
    p.add_argument('--outdir', default='figures/fig6')
    p.add_argument('--stem', default='fig6C_projection')
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

    conf = pd.read_csv(indir / 'confusion_coarse_celltype.tsv', sep='\t', index_col=0)
    land = pd.read_csv(indir / 'query_ne_landing.tsv', sep='\t')

    conf = conf.reindex(index=[r for r in ROW_ORDER if r in conf.index],
                        columns=[c for c in COL_ORDER if c in conf.columns]).fillna(0)
    frac = conf.div(conf.sum(1), axis=0)

    from _utils import setup_style, savefig
    plt = setup_style()
    fig, (axC, axL) = plt.subplots(1, 2, figsize=figsize(8.6, 3.8),
                                   gridspec_kw={'width_ratios': [1.15, 1]})

    # ---- left: confusion heatmap -------------------------------------------
    im = axC.imshow(frac.values, cmap='Blues', vmin=0, vmax=1, aspect='auto')
    axC.set_xticks(range(frac.shape[1]))
    axC.set_xticklabels([PRETTY.get(c, c) for c in frac.columns], rotation=45,
                        ha='right', fontsize=7.5)
    axC.set_yticks(range(frac.shape[0]))
    axC.set_yticklabels([PRETTY.get(r, r) for r in frac.index], fontsize=7.5)
    axC.set_xlabel('transferred LTL331 label', fontsize=8.5)
    axC.set_ylabel('clinical original lineage', fontsize=7.5)
    axC.set_title('projection label transfer', fontsize=8)
    for i in range(frac.shape[0]):
        for j in range(frac.shape[1]):
            v = frac.values[i, j]
            if v > 0.01:
                axC.text(j, i, f'{v:.0%}', ha='center', va='center', fontsize=7.5,
                         color='white' if v > 0.55 else 'black')
    cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.02)
    cb.set_label('row fraction', fontsize=8.5); cb.ax.tick_params(labelsize=7.5)

    # ---- right: per-cohort NE landing --------------------------------------
    piv = land.pivot_table(index='_ds', columns='_t', values='n', aggfunc='sum').fillna(0)
    piv = piv.reindex(columns=[c for c in COL_ORDER if c in piv.columns]).fillna(0)
    pf = piv.div(piv.sum(1), axis=0)
    cohorts = [c for c in ['Gao', 'Li'] if c in pf.index] + \
              [c for c in pf.index if c not in ('Gao', 'Li')]
    pf = pf.reindex(cohorts)
    bottom = np.zeros(len(pf))
    x = np.arange(len(pf))
    for col in pf.columns:
        axL.bar(x, pf[col].values, bottom=bottom, width=0.6,
                color=TCOLORS.get(col, '#999999'), label=PRETTY.get(col, col))
        bottom += pf[col].values
    axL.set_xticks(x); axL.set_xticklabels(pf.index, fontsize=7.5)
    axL.set_ylabel('fraction of clinical NE cells', fontsize=7.5)
    axL.set_ylim(0, 1)
    axL.set_title('where clinical NE cells land', fontsize=8)
    axL.legend(frameon=False, fontsize=7.5, ncol=1, loc='center left',
               bbox_to_anchor=(1.0, 0.5))
    axL.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    savefig(fig, outdir, args.stem, ('pdf', 'png'))
    plt.close(fig)
    ne_row = frac.loc['NE'] if 'NE' in frac.index else None
    if ne_row is not None:
        ne_to_ne = ne_row.get('NE_ASCL1pos', 0) + ne_row.get('NE_ASCL1neg', 0)
        print(f'[6C] PRAD->PRAD={frac.loc["PRAD","PRAD"]:.0%} | NE->NE={ne_to_ne:.0%}')


if __name__ == '__main__':
    main()
