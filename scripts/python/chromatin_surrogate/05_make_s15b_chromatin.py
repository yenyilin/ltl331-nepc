#!/usr/bin/env python3
"""
05_make_s15b_chromatin.py — assemble the Supplementary Fig. S15B chromatin-surrogate panel.

Parses every HOMER knownResults.txt under $OUT/motif/*/, maps HOMER motif names to
our regulon TF set, and renders a heatmap of motif enrichment (-log10 p) for our TFs
(rows) across the NEPC-vs-adeno contrasts (columns). This is the figure that answers
a natural question: "without matched ATAC, do our nominated regulators show concordant
motif enrichment in clinical/PDX NEPC chromatin?"

HOMER knownResults.txt columns (tab, with header):
  Motif Name | Consensus | P-value | Log P-value(ln) | q-value | #Target.. | %Target | #Bg.. | %Bg
We use Log P-value (natural log, negative); -log10 p = -lnP / ln(10).

Run:  python 05_make_s15b_chromatin.py            (reads $OUT, $REGULON_TFS)
      python 05_make_s15b_chromatin.py --motif-root <dir> --outdir figures
"""
import argparse, os, re, math
from pathlib import Path
import numpy as np
import pandas as pd

# regulon TF -> regex matched against the HOMER "Motif Name" field.
# Ordered as two blocks: the NE program (top) and the neural-crest / EMT tier
# (bottom). The block split is drawn as a divider in the panel so the large
# strength gap between the NE headliners and the weaker NC-specifier motifs is
# transparent rather than hidden by the color cap.
TF_PATTERNS = {
    # --- NE program ---
    'MSX1':    r'\bMSX1?\b',
    'ASCL1':   r'ASCL1',
    'NEUROD1': r'NEUROD1',
    'FOXA2':   r'FOXA2',
    'FOXA1':   r'FOXA1',
    'NKX2-1':  r'NKX2[._-]?1',
    'POU3F2':  r'POU3F2|BRN2',
    'INSM1':   r'INSM1',
    'SOX2':    r'\bSOX2\b',
    'ONECUT2': r'ONECUT2?',
    'REST':    r'\bREST\b|NRSF',
    'RUNX3':   r'RUNX3',
    # --- neural-crest / EMT tier ---
    'SNAI2':   r'SNAI2|SLUG',
    'SNAI1':   r'SNAI1|SNAIL1',
    'TWIST1':  r'TWIST1?',
    'ZEB1':    r'\bZEB1\b',
    'ZEB2':    r'\bZEB2\b',
    'FOXD3':   r'FOXD3',
    'SOX10':   r'SOX10',
    'TFAP2C':  r'AP-?2gamma|TFAP2C',
    'PAX7':    r'PAX7',
}
# first TF of the neural-crest / EMT tier (used to draw the block divider)
NC_TIER_START = 'SNAI2'
LN10 = math.log(10)


def parse_known(path):
    """Return dict TF -> best -log10(p) found in this knownResults.txt."""
    df = pd.read_csv(path, sep='\t')
    # the log-p column name varies slightly across HOMER versions
    logcol = next((c for c in df.columns if re.search(r'log\s*p[- ]?value', c, re.I)), None)
    namecol = df.columns[0]
    best = {}
    for tf, pat in TF_PATTERNS.items():
        hit = df[df[namecol].astype(str).str.contains(pat, case=False, regex=True)]
        if hit.empty:
            continue
        lnp = hit[logcol].astype(float).min()   # most significant (most negative ln p)
        best[tf] = -lnp / LN10                   # -log10 p
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--motif-root', default=os.environ.get('OUT', 'out') + '/motif')
    ap.add_argument('--tfs', nargs='*',
                    default=os.environ.get('REGULON_TFS', '').split() or list(TF_PATTERNS))
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='s15b_chromatin_motif')
    ap.add_argument('--vmax', type=float, default=30.0, help='cap on -log10 p for color')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author panel at this final width in mm (S15 cell target) so it scales ~1:1')
    ap.add_argument('--figheight-mm', type=float, default=None)
    args = ap.parse_args()

    runs = sorted(Path(args.motif_root).glob('*/knownResults.txt'))
    if not runs:
        raise SystemExit(f'no knownResults.txt under {args.motif_root} — run 03 first')

    data = {}   # contrast label -> {TF: -log10p}
    for kr in runs:
        label = kr.parent.name
        data[label] = parse_known(kr)
        print(f'[parsed] {label}: {len(data[label])} of our TFs matched')

    tfs = [t for t in args.tfs if t in TF_PATTERNS]
    cols = list(data.keys())
    M = pd.DataFrame(index=tfs, columns=cols, dtype=float)
    for c in cols:
        for t in tfs:
            M.loc[t, c] = data[c].get(t, np.nan)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    M.to_csv(outdir / f'{args.stem}_data.tsv', sep='\t')
    print(f'[data] wrote {args.stem}_data.tsv\n{M.round(1).to_string()}')

    import matplotlib
    matplotlib.use('Agg'); matplotlib.rcParams['pdf.fonttype'] = 42
    import matplotlib.pyplot as plt
    D = M.astype(float).clip(upper=args.vmax)

    def figsize(dw, dh):
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh: return (fw / 25.4, fh / 25.4)
        if fw:        return (fw / 25.4, dh * (fw / 25.4) / dw)
        if fh:        return (dw * (fh / 25.4) / dh, fh / 25.4)
        return (dw, dh)

    fig, ax = plt.subplots(figsize=figsize(0.7 * len(cols) + 3.5, 0.32 * len(tfs) + 1.5))
    im = ax.imshow(D.values, aspect='auto', cmap='Reds', vmin=0, vmax=args.vmax)
    short = {'NECRE_vs_AdCRE': 'NE vs Ad\nCRE', 'baca_FOXA1_NEPC_vs_adeno': 'FOXA1\nNE vs Ad'}
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([short.get(c, c) for c in cols], rotation=0, ha='center', fontsize=7)
    ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs, fontsize=7.5)
    # annotate each cell with the TRUE (unclipped) -log10 p, so the color cap does
    # not hide the >vmax headliner values or the weaker NC-tier gradient
    Mv = M.astype(float).values
    for i in range(len(tfs)):
        for j in range(len(cols)):
            vtrue = Mv[i, j]
            if not np.isnan(vtrue):
                ax.text(j, i, f'{vtrue:.0f}', ha='center', va='center', fontsize=7,
                        color='white' if D.values[i, j] > args.vmax * 0.55 else 'black')
    # divider between the NE program and the neural-crest / EMT tier
    if NC_TIER_START in tfs:
        b = tfs.index(NC_TIER_START)
        ax.axhline(b - 0.5, color='0.25', lw=1.2)   # NE program (top) | NC/EMT tier (bottom)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, aspect=12, pad=0.02)
    cb.set_label('−log10 p (HOMER)', fontsize=7)
    cb.ax.tick_params(labelsize=7)
    ax.set_title('Motif enrichment,\nNEPC vs adeno CREs', fontsize=8.5)
    ax.tick_params(length=0)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}')


if __name__ == '__main__':
    main()
