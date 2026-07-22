#!/usr/bin/env python3
"""
cnv_driver_locus_tracks.py — copy-number state of the canonical prostate/NEPC
driver loci across the castration timecourse, in BOTH modalities (scRNA inferCNV
+ bulk WES). Supports two linked arguments:

  (1) Permissive-genome thesis: the NEPC-enabling lesions (chr8q/MYC gain,
      RB1 / TP53 / PTEN losses, MYCN) are present ALREADY at pre-castration and
      INVARIANT through the PRAD->NEPC transition -> transdifferentiation is
      epigenetic/transcriptional on a fixed genome, not driven by new CNAs.
  (2) The AR point: the AR locus (chrX) is copy-number STABLE across castration
      even though AR *expression* collapses -> AR loss is transcriptional
      silencing, not genomic deletion. (Adds mechanistic depth to the trajectory
      model, complementing the WES-inferCNV copy-number concordance.)

Reads the two tidy calls TSVs (nexus_cnv_loader / infercnv_loader schema) and,
for each driver gene, reports the copy state of any segment overlapping the gene
body (max-|value| if several; neutral 0 if none). Emits a gene x timepoint table
per modality + a paired heatmap.

Example:
  python cnv_driver_locus_tracks.py \\
      --infercnv-calls data/infercnv_calls.tsv \\
      --wes-calls data/nexus_calls.tsv \\
      --outdir figures/cnv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Driver loci, hg19 (GRCh37) gene bodies; dir = expected NEPC-associated change.
DRIVERS = [
    ('MYC',  'chr8',  128748315, 128753680, 'gain'),   # 8q24 amplicon
    ('MYCN', 'chr2',   15940550,  15947007, 'gain'),   # 2p24 (NEPC)
    ('AR',   'chrX',   66764465,  66950461, 'AR'),      # castration target
    ('RB1',  'chr13',  48877883,  49056122, 'loss'),    # NEPC tumor suppressor
    ('TP53', 'chr17',   7565097,   7590856, 'loss'),    # NEPC tumor suppressor
    ('PTEN', 'chr10',  89623195,  89728532, 'loss'),
    ('CDH1', 'chr16',  68771128,  68869444, 'loss'),    # 16q22 (the chr16q loss in S3)
]

# scRNA inferCNV timepoints (samples.groups) — already named by timepoint.
IC_TP_ORDER = ['preCx', '2wk', '4wk', '8wk', '12wk', '16wk', '20wk', '22wk']

# WES/Nexus sample -> timepoint. Keys use the public sample names from
# nexus_cnv_loader.py (no normalization-reference suffix).
# ONLY the four SAME-STUDY matched WES timepoints (preCx, 16, 20, 22 wk) are used:
# each was newly generated in this study and normalized to its matched germline
# normal. The 8- and 12-wk tracks are LTL331-5 aCGH from a SEPARATE xenograft
# generation and a different platform/normalization reference — they are NOT
# matched WES and are excluded here (consistent with the S6 concordance footnote).
WES_TP = {
    'LTL331_preCX1': 'preCx',
    'LTL331_16wk1':  '16wk',
    'LTL331_20wk':   '20wk',
    'LTL331_22wk2':  '22wk',
}
WES_TP_ORDER = ['preCx', '16wk', '20wk', '22wk']


def locus_state(calls, sample, chrom, start, end):
    """Copy state of the segment(s) overlapping [start,end] on chrom for sample;
    max-|value| if several, 0 (neutral) if none."""
    sub = calls[(calls['sample'].astype(str) == sample) &
                (calls['Chromosome'] == chrom) &
                (calls['Start'] < end) & (calls['End'] > start)]
    if sub.empty:
        return 0
    v = sub['Value'].to_numpy()
    return int(v[np.argmax(np.abs(v))])


def build_table(calls, sample_to_tp, tp_order):
    rows = {}
    inv = {tp: s for s, tp in sample_to_tp.items()}
    for gene, chrom, s, e, _ in DRIVERS:
        rows[gene] = [locus_state(calls, inv[tp], chrom, s, e) if tp in inv else np.nan
                      for tp in tp_order]
    return pd.DataFrame(rows, index=tp_order).T  # genes x timepoints


def draw(ax, mat, title):
    A = mat.to_numpy(dtype=float)
    im = ax.imshow(A, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, fontsize=7, rotation=45, ha='right')
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_title(title, fontsize=8)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if np.isnan(A[i, j]):
                ax.text(j, i, '·', ha='center', va='center', fontsize=8, color='0.5')
            else:
                ax.text(j, i, f'{int(A[i, j]):+d}'.replace('+0', '0'),
                        ha='center', va='center', fontsize=6.5,
                        color='white' if abs(A[i, j]) >= 2 else 'black')
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--infercnv-calls', default='data/infercnv_calls.tsv')
    ap.add_argument('--wes-calls', default='data/nexus_calls.tsv')
    ap.add_argument('--outdir', default='figures/cnv')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    ic = pd.read_csv(args.infercnv_calls, sep='\t', dtype={'sample': str})
    we = pd.read_csv(args.wes_calls, sep='\t', dtype={'sample': str})

    ic_tab = build_table(ic, {tp: tp for tp in IC_TP_ORDER}, IC_TP_ORDER)
    we_tab = build_table(we, WES_TP, WES_TP_ORDER)
    ic_tab.to_csv(outdir / 'driver_locus_inferCNV.tsv', sep='\t')
    we_tab.to_csv(outdir / 'driver_locus_WES.tsv', sep='\t')

    # invariance check: per gene, is the DIRECTION (sign: loss/neutral/gain)
    # constant across timepoints? Amplitude wobble (1 vs 2 copies) is not a
    # biological state change; a sign flip would be. Reports the preCx state too.
    def sign_inv(tab):
        out = {}
        for g in tab.index:
            row = tab.loc[g].dropna()
            signs = set(np.sign(row.to_numpy()).astype(int))
            pre = tab.loc[g].get('preCx', np.nan)
            out[g] = (len(signs) == 1, {-1: 'loss', 0: 'neutral', 1: 'gain'}.get(
                int(np.sign(pre)) if not np.isnan(pre) else 9, '?'))
        return out
    ici, wei = sign_inv(ic_tab), sign_inv(we_tab)
    print('[driver-locus direction-invariance]  (preCx state + sign-constant => preset & permissive)')
    for g in ic_tab.index:
        print(f'  {g:5s}  WES preCx={wei[g][1]:7s} {"sign-CONST" if wei[g][0] else "sign-VARIES":11s}'
              f' | inferCNV {"sign-CONST" if ici[g][0] else "sign-VARIES":11s}')

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.2),
                             gridspec_kw={'width_ratios': [len(IC_TP_ORDER), len(WES_TP_ORDER)]})
    draw(axes[0], ic_tab, 'scRNA inferCNV (per timepoint)')
    im = draw(axes[1], we_tab, 'bulk WES (matched: preCx, 16, 20, 22 wk)')
    fig.suptitle('Driver-locus copy number across the castration timecourse — '
                 'preset at pre-Cx and invariant (AR locus stable despite AR-expression loss)',
                 fontsize=9)
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02, ticks=[-2, -1, 0, 1, 2])
    cb.set_label('copy state', fontsize=7); cb.ax.tick_params(labelsize=6)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'driver_locus_tracks.{fmt}', dpi=300, bbox_inches='tight')
    print(f'[ok] wrote {outdir}/driver_locus_tracks.{{{",".join(args.plot_format)}}} + 2 TSVs')


if __name__ == '__main__':
    main()
