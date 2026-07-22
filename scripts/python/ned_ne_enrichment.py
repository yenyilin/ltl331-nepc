#!/usr/bin/env python3
"""
ned_ne_enrichment.py — Fig 6D: neuroendocrine epithelium is enriched in patients WITH
clinical neuroendocrine differentiation (NED) vs without, in both clinical cohorts.

Computes the NE-vs-rest x NED-vs-nonNED 2x2 Fisher exact odds ratio + p per clinical
cohort straight from the coarse-label-by-NED counts (no object needed), writes the stats
TSV the Fig 6 layout expects (validation_ne_enrichment.tsv), and renders a forest plot.

Woolf 95% CI on the log-odds (0.5 continuity if any zero cell).

Example:
  python scripts/ned_ne_enrichment.py \
      --counts data/integration_qc/validation_coarse_by_ned.tsv \
      --outdir figures/fig6 --stats-out data/integration_qc/validation_ne_enrichment.tsv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def woolf_ci(a, b, c, d):
    """95% CI for OR via Woolf log method; 0.5 continuity if any zero."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    logor = np.log((a * d) / (b * c))
    return float(np.exp(logor - 1.96 * se)), float(np.exp(logor + 1.96 * se))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--counts', default='data/integration_qc/validation_coarse_by_ned.tsv')
    p.add_argument('--ne-label', default='NE')
    p.add_argument('--outdir', default='figures/fig6')
    p.add_argument('--stem', default='fig6D_ned_enrichment')
    p.add_argument('--stats-out', default='data/integration_qc/validation_ne_enrichment.tsv')
    p.add_argument('--figwidth-mm', type=float, default=None,
                   help='author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa')
    p.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm')
    args = p.parse_args()

    outdir = Path(args.outdir)
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

    df = pd.read_csv(args.counts, sep='\t')
    # clinical cohorts only (have NED / non-NED); LTL331 is PDX (single status) -> skip
    clinical = df[df['ned_status'].isin(['NED', 'non-NED'])]
    rows = []
    for cohort, sub in clinical.groupby('dataset'):
        piv = sub.pivot_table(index='coarse_celltype', columns='ned_status', values='n',
                              aggfunc='sum').fillna(0)
        if not {'NED', 'non-NED'}.issubset(piv.columns) or args.ne_label not in piv.index:
            continue
        a = piv.loc[args.ne_label, 'NED']                    # NE & NED
        c = piv.loc[args.ne_label, 'non-NED']                # NE & non-NED
        b = piv.drop(args.ne_label)['NED'].sum()             # non-NE & NED
        d = piv.drop(args.ne_label)['non-NED'].sum()         # non-NE & non-NED
        orr, pval = fisher_exact([[a, b], [c, d]])
        lo, hi = woolf_ci(a, b, c, d)
        rows.append({'cohort': cohort, 'NE_NED': int(a), 'nonNE_NED': int(b),
                     'NE_nonNED': int(c), 'nonNE_nonNED': int(d),
                     'odds_ratio': round(orr, 2), 'ci_low': round(lo, 2),
                     'ci_high': round(hi, 2), 'p_value': pval})
    stats = pd.DataFrame(rows)
    Path(args.stats_out).parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.stats_out, sep='\t', index=False)
    print('=== NED -> NE enrichment ===')
    print(stats.to_string(index=False))

    # ---- forest plot --------------------------------------------------------
    from _utils import setup_style, savefig
    plt = setup_style()
    fig, ax = plt.subplots(figsize=figsize(4.3, 1.0 + 0.5 * len(stats)))
    y = np.arange(len(stats))[::-1]
    ax.errorbar(stats['odds_ratio'], y,
                xerr=[stats['odds_ratio'] - stats['ci_low'],
                      stats['ci_high'] - stats['odds_ratio']],
                fmt='o', color='#b2182b', ecolor='0.4', elinewidth=1, capsize=3,
                markersize=6)
    ax.axvline(1, color='0.5', lw=0.7, ls=':')
    ax.set_xscale('log')
    for yi, r in zip(y, stats.itertuples()):
        ax.text(r.ci_high * 1.15, yi, f'OR {r.odds_ratio:g} ({r.ci_low:g}–{r.ci_high:g})',
                va='center', fontsize=8.5)
    ax.set_yticks(y); ax.set_yticklabels(stats['cohort'], fontsize=9)
    ax.set_ylim(-0.6, len(stats) - 0.4)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_xlabel('odds ratio — NE epithelium in NED vs non-NED (log scale)', fontsize=8.5)
    ax.set_title('clinical NE differentiation ↔ NE state', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    savefig(fig, outdir, args.stem, ('pdf', 'png'))
    plt.close(fig)


if __name__ == '__main__':
    main()
