#!/usr/bin/env python3
"""
plot_scenic_regulon_map.py — Supplementary Fig. S13 (RSS): SCENIC regulon *specificity*
(how exclusively each regulon marks a cluster) across the adenocarcinoma -> EMT ->
neuroendocrine clusters. This is the SPECIFICITY map; the regulon-ACTIVITY map (Fig 5A,
per-cluster mean AUCell) is plot_scenic_auc_map.py — the two orderings agree for
cluster-restricted regulators and diverge for broadly active ones.

Input is the per-cluster RSS *rank* table (data/scenic/all_all_rss_rank.tsv: rows =
regulons like ASCL1_(+), cols = clusters, value = RSS rank where 1 = most cluster-specific).
Rank is converted to a specificity score in [0,1] (1 = rank 1) for the heatmap. Rows are
ordered diagonally by the cluster where each regulon is most specific, so PRAD regulons sit
top-left and NE regulons (ASCL1 etc.) bottom-right along the canonical axis.

Shows the union of the top-N most cluster-specific regulons per cluster, plus any --highlight
TFs present (bolded). SCENIC regulons are activating _(+) or repressing _(-); BOTH are shown
by default and sign-labeled (a TF's (+) and (-) regulons often peak in different clusters,
e.g. FOXA2 (+) in NE vs FOXA2 (-) in AR-low). Use --exclude-repressors for activating-only.
Clusters follow data/cluster_order.json.

Example:
  python plot_scenic_regulon_map.py --rss data/scenic/all_all_rss_rank.tsv \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --na-tfs NKX2-1 ONECUT2 POU3F2 INSM1 \\
      --top 5 --outdir figures --stem s13_scenic_rss   # top-5 per cluster, matching Fig 5A
"""
import argparse, json, re
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_HL = ['MSX1', 'ASCL1', 'FOXA2', 'NEUROD1', 'INSM1', 'TP63', 'SNAI2',
              'RFX2', 'RFX3', 'AR', 'FOXA1', 'GATA2', 'HOXB13']


def base_tf(name):
    """TF name without the +/- sign (for matching/highlight)."""
    return re.sub(r'_\([+-]\)$', '', str(name))


def sign_of(name):
    return '-' if str(name).endswith('_(-)') else '+'


def disp(name, show_sign):
    """display label; keep the sign only when repressors are also shown."""
    b = base_tf(name)
    return f'{b} ({sign_of(name)})' if show_sign else b


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rss', default='data/scenic/all_all_rss_rank.tsv')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--cluster-colors', default=None)
    ap.add_argument('--top', type=int, default=5, help='top-N specific regulons per cluster')
    ap.add_argument('--regulons', nargs='*', default=None,
                    help='explicit TF list (overrides --top selection)')
    ap.add_argument('--highlight', nargs='+', default=DEFAULT_HL)
    ap.add_argument('--na-tfs', nargs='*', default=None,
                    help='canonical TFs to append as greyed "N.A. in SCENIC" rows at the '
                         'bottom (only those with NO SCENIC regulon are drawn). Visual '
                         'argument for Fig 5B/decoupleR.')
    ap.add_argument('--exclude-repressors', action='store_true',
                    help='show only activating _(+) regulons (default: show BOTH _(+) and _(-), '
                         'sign-labeled — they peak in different clusters)')
    ap.add_argument('--cmap', default='magma_r')
    ap.add_argument('--figwidth', type=float, default=None)
    ap.add_argument('--figheight', type=float, default=None)
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig5A_scenic_regulons')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    rss = pd.read_csv(args.rss, sep='\t', index_col=0)
    rss.columns = rss.columns.map(str)
    n_repr = int(sum(str(r).endswith('_(-)') for r in rss.index))
    if args.exclude_repressors:
        rss = rss.loc[[r for r in rss.index if not str(r).endswith('_(-)')]]
        print(f'[sign] activating (+) regulons only; excluded {n_repr} repressing (-) regulons')
    else:
        print(f'[sign] showing BOTH (+) and (-) regulons, sign-labeled ({n_repr} repressors '
              f'included; use --exclude-repressors for activating-only)')
    N = rss.shape[0]
    cats = list(rss.columns)
    if Path(args.cluster_order).exists():
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        cols = [c for c in co if c in cats] + [c for c in cats if c not in co]
    else:
        cols = sorted(cats, key=lambda s: int(s) if s.isdigit() else 1e9)

    hl_set = {h.upper() for h in args.highlight}
    sel = set()
    if args.regulons:
        want = {r.upper() for r in args.regulons}
        sel = {r for r in rss.index if base_tf(r).upper() in want}
    else:
        for c in cols:
            sel.update(rss[c].nsmallest(args.top).index)
    # always add present highlight regulons
    hl_rows = {r for r in rss.index if base_tf(r).upper() in hl_set}
    sel |= hl_rows
    sel = list(sel)
    if not sel:
        raise SystemExit('no regulons selected')
    missing_hl = sorted(hl_set - {base_tf(r).upper() for r in hl_rows})
    if missing_hl:
        print(f'[note] highlight TFs with no SCENIC regulon (show these in Fig 5B/decoupleR): '
              f'{missing_hl}')

    M = rss.loc[sel, cols]
    S = 1.0 - (M - 1) / (N - 1)                      # specificity in [0,1]
    pos = {c: i for i, c in enumerate(cols)}
    best_c = M.idxmin(axis=1); best_r = M.min(axis=1)
    rowkey = sorted(sel, key=lambda r: (pos.get(best_c[r], 99), best_r[r]))
    S = S.loc[rowkey]

    # optional greyed "N.A. in SCENIC" block (canonical TFs absent from the regulon set)
    na_label = {}
    if args.na_tfs:
        present = {base_tf(r).upper() for r in S.index}
        na_draw = [t for t in args.na_tfs if t.upper() not in present]
        na_skip = [t for t in args.na_tfs if t.upper() in present]
        if na_skip:
            print(f'[na] --na-tfs that ARE SCENIC regulons (not greyed, shown above): {na_skip}')
        if na_draw:
            idx = [f'__NA__{t.upper()}' for t in na_draw]
            na_label = {i: f'{t} (N.A.)' for i, t in zip(idx, na_draw)}
            S = pd.concat([S, pd.DataFrame(np.nan, index=idx, columns=S.columns)])
            print(f'[na] greyed "N.A. in SCENIC" rows appended: {na_draw}')

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt

    nrow, ncol = S.shape
    W = args.figwidth or (0.34 * ncol + 3.0)
    H = args.figheight or (0.17 * nrow + 1.8)
    fig, ax = plt.subplots(figsize=(W, H))
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad('0.8')                              # NaN (N.A. in SCENIC) -> grey
    im = ax.imshow(S.values, aspect='auto', cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(cols, fontsize=9, rotation=90, va='top')
    rowkey_full = list(S.index)                      # rowkey + appended N.A. rows
    ax.set_yticks(range(nrow))
    show_sign = not args.exclude_repressors
    ylabs = [na_label.get(r, disp(r, show_sign)) for r in rowkey_full]
    ax.set_yticklabels(ylabs, fontsize=6.5)
    for t, r in zip(ax.get_yticklabels(), rowkey_full):
        if r in na_label:
            t.set_color('0.45'); t.set_style('italic')
        elif base_tf(r).upper() in hl_set:
            t.set_fontweight('bold')
    if na_label:
        ax.axhline(nrow - len(na_label) - 0.5, color='0.6', lw=0.6, ls='--')
    colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if (args.cluster_colors and Path(args.cluster_colors).exists()) else {}
    if colors:
        for t, c in zip(ax.get_xticklabels(), cols):
            t.set_color(colors.get(c, 'black'))
    ax.set_xlabel('cluster (adenocarcinoma → neuroendocrine)', fontsize=9)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.4, aspect=28, pad=0.02)
    cb.set_label('regulon specificity (RSS)', fontsize=10); cb.ax.tick_params(labelsize=8)
    ax.set_title('SCENIC regulon specificity across NEtD', fontsize=11)
    if na_label:
        ax.text(0.0, -0.14, 'N.A. = not recovered by SCENIC (absent from the coexpression/specificity step)',
                transform=ax.transAxes, fontsize=6, color='0.45', style='italic', va='top')

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}  '
          f'({nrow} regulons x {ncol} clusters)')


if __name__ == '__main__':
    main()
