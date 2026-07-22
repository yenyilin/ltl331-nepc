#!/usr/bin/env python3
"""
plot_shrestha_chromatin.py — main Fig 5D chromatin anchor and its full supplementary
companion (S15A).

Clinical-tissue chromatin validation: regulon TF footprint occupancy across mCRPC
transcriptional subtypes from Shrestha et al. 2024 Cancer Research (WCDT, n=70
clinical ATAC-seq biopsies), Supplementary Table S5. Shows that our nominated
regulators (ASCL1, NEUROD1, POU3F2, the novel MSX1, plus NHLH1/cl9 and the cl17
EMT TFs) have AR-NE+ (NEPC)-specific chromatin accessibility in patient tissue.

--curated  -> MAIN Fig 5D: the story rows only (NE drivers + EMT/bridge), stem fig5D.
default    -> full 17-TF atlas (adds AR-axis + bHLH-partner context) for supp S15A.

Input: data/chromatin_surrogate/shrestha2024_S5_regulon_occupancy.tsv
Open supplementary data — no controlled-access (EGA) download required.

Honesty rendered IN the figure: per-subtype n in the column labels; TFs absent from
the 203-TF footprint set (NKX2-1/FOXA2/SOX2/REST/RUNX3) listed as a footnote, not
silently dropped.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# subtype column order: adenocarcinoma -> NE; with per-subtype n (Supp Table S1)
SUBTYPE_ORDER = ['AR+NE-', 'ARlowNE-', 'AR+NE+', 'AR-NE-', 'AR-NE+']
SUBTYPE_N = {'AR+NE-': 26, 'ARlowNE-': 32, 'AR+NE+': 2, 'AR-NE-': 6, 'AR-NE+': 4}

# TF row groups (PRAD/broad -> bHLH -> EMT -> NE drivers; NE drivers last/bottom)
GROUPS = [
    ('AR-axis / broad', ['FOXA1', 'ONECUT2']),
    ('EMT / bridge',    ['SNAI2', 'ZEB1', 'TWIST1']),
    ('bHLH partners',   ['TCF3', 'TCF4', 'ID4']),
    ('NE drivers',      ['ASCL1', 'NEUROD1', 'POU3F2', 'MSX1', 'NHLH1', 'POU4F3', 'INSM1',
                         'NKX2-1', 'FOXA2', 'SOX2', 'REST', 'RUNX3']),
]
HL = {'ASCL1', 'NEUROD1', 'POU3F2', 'MSX1', 'NHLH1'}   # bolded


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--occupancy',
                    default='data/chromatin_surrogate/shrestha2024_S5_regulon_occupancy.tsv')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig5C_chromatin_shrestha')
    ap.add_argument('--curated', action='store_true',
                    help='MAIN Fig 5D view: keep only the story TF groups '
                         '(NE drivers + EMT/bridge). Default (no flag) = full 17-TF '
                         'atlas for the supplementary companion (S15A).')
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

    groups_use = ([g for g in GROUPS if g[0] in ('NE drivers', 'EMT / bridge')]
                  if args.curated else GROUPS)
    if args.curated:
        print('[curated] MAIN 5D view: NE drivers + EMT/bridge only '
              '(full atlas is the default, used for S15A)')

    df = pd.read_csv(args.occupancy, sep='\t')
    df['tf_name'] = df['tf_name'].astype(str)
    df = df.set_index('tf_name')

    rows, groups, absent = [], {}, []
    for gname, tfs in groups_use:
        for t in tfs:
            if t in df.index and str(df.loc[t, 'AR_NE_pos_rank_of203']) != 'absent':
                rows.append(t); groups[t] = gname
            else:
                absent.append(t)
    M = df.loc[rows, SUBTYPE_ORDER].astype(float)

    import matplotlib
    matplotlib.use('Agg'); matplotlib.rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt

    nrow, ncol = M.shape
    fig, ax = plt.subplots(figsize=figsize(0.95 * ncol + 3.0, 0.34 * nrow + 2.2))
    vmax = float(np.nanpercentile(M.values, 98)) or 1.0
    im = ax.imshow(M.values, aspect='auto', cmap='Reds', vmin=0, vmax=vmax)

    collabels = [f'{s} (n={SUBTYPE_N[s]})' for s in SUBTYPE_ORDER]
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(collabels, fontsize=7, rotation=40, ha='right',
                       rotation_mode='anchor')
    # mark the NEPC column
    ne_i = SUBTYPE_ORDER.index('AR-NE+')
    ax.get_xticklabels()[ne_i].set_color('#b30000'); ax.get_xticklabels()[ne_i].set_fontweight('bold')
    ax.set_yticks(range(nrow)); ax.set_yticklabels(rows, fontsize=8)
    for t, r in zip(ax.get_yticklabels(), rows):
        if r in HL: t.set_fontweight('bold')

    # value annotations
    for i in range(nrow):
        for j in range(ncol):
            v = M.values[i, j]
            if v > 0:
                ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=6.5,
                        color='white' if v > vmax * 0.55 else '0.2')
    # group separators
    seen, bnd = [], []
    for i, r in enumerate(rows):
        if groups[r] not in seen:
            seen.append(groups[r])
            if i > 0: bnd.append(i)
    for b in bnd: ax.axhline(b - 0.5, color='0.5', lw=0.8)
    tr = ax.get_yaxis_transform(); starts = [0] + bnd + [nrow]
    for gi in range(len(seen)):
        mid = (starts[gi] + starts[gi+1] - 1) / 2
        # wrap at '/' or space so the vertical label fits a short group span (2-3 rows)
        lab = seen[gi].replace(' / ', '\n').replace(' ', '\n')
        ax.text(1.02, mid, lab, transform=tr, ha='left', va='center',
                fontsize=6.5, rotation=90, color='0.3', linespacing=0.9)

    cb = fig.colorbar(im, ax=ax, shrink=0.45, aspect=12, pad=0.16)
    cb.set_label('TF footprint occupancy score (Shrestha 2024, Supp S5)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    # keep each title line short: a long one-liner floors the panel width under bbox='tight'
    ax.set_title('Clinical mCRPC chromatin (WCDT ATAC, n=70):\n'
                 'regulon TF footprint is NEPC-specific', fontsize=9)
    ax.tick_params(length=0)
    foot = ('AR-NE+ subtype n=4 (descriptive). Absent from the 203 differential-footprint TFs: '
            + ', '.join(absent) + '.')
    fig.text(0.0, -0.13, foot, fontsize=6.5, color='0.35', wrap=True)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}  ({nrow} TFs x {ncol} subtypes)')
    print(f'     absent from footprint set (footnoted): {absent}')


if __name__ == '__main__':
    main()
