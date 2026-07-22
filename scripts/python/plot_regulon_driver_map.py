#!/usr/bin/env python3
"""
plot_regulon_driver_map.py — Fig 5B: comparative TF-activity map (decoupleR / CollecTRI)
that places MSX1 alongside the canonical NEPC drivers (NKX2-1, ONECUT2, POU3F2/BRN2) plus
other lineage TFs.

decoupleR (run via scripts/decoupler_tf_activity.py) scores TF activity per cell from a
curated TF->target prior, bypassing SCENIC's coexpression+motif filter, so TFs SCENIC
pruned (NKX2-1, ONECUT2, POU3F2 — see review/scenic_rss_findings.md) are recovered here.

Input: tf_activity_per_cluster.tsv from decoupler_tf_activity.py (rows = TF, cols = cluster,
value = ULM/MLM activity estimate). Rows shown in a curated PRAD -> EMT -> NE order; color =
per-TF z-score across clusters (diverging) by default; canonical NEPC drivers + MSX1 bolded.

Example:
  python plot_regulon_driver_map.py \\
      --activity data/decoupler_tf/tf_activity_per_cluster.tsv \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --outdir figures --stem fig5B_tf_drivers
"""
import argparse, json
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

# curated TFs grouped by program (order = PRAD -> EMT -> NE); only those present are shown.
DEFAULT_TFS = OrderedDict([
    ('AR-axis / PRAD', ['AR', 'FOXA1', 'NKX3-1', 'HOXB13', 'GATA2']),
    ('EMT',            ['SNAI2', 'SNAI1', 'TWIST1', 'ZEB1', 'TEAD1']),
    ('NE drivers',     ['ASCL1', 'NEUROD1', 'INSM1', 'FOXA2', 'NKX2-1', 'ONECUT2',
                        'POU3F2', 'SOX2', 'MSX1']),
    ('NE repressor',   ['REST']),
])
# canonical NEPC drivers + the novel nomination -> bolded on the y-axis
HL = {'NKX2-1', 'ONECUT2', 'POU3F2', 'ASCL1', 'NEUROD1', 'INSM1', 'FOXA2', 'MSX1'}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--activity', default='data/decoupler_tf/tf_activity_per_cluster.tsv')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--cluster-colors', default=None)
    ap.add_argument('--tfs', nargs='*', default=None,
                    help='explicit TF order (overrides the curated grouped default)')
    ap.add_argument('--scale', choices=['zscore', 'none'], default='zscore')
    ap.add_argument('--vlim', type=float, default=2.0)
    ap.add_argument('--highlight', nargs='+', default=sorted(HL))
    ap.add_argument('--figwidth', type=float, default=None)
    ap.add_argument('--figheight', type=float, default=None)
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig5B_tf_drivers')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    act = pd.read_csv(args.activity, sep='\t', index_col=0)
    act.columns = act.columns.map(str)
    act.index = act.index.map(str)
    # decoupler may write clusters as rows; orient so rows = TF, cols = cluster
    cats_json = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])] \
        if Path(args.cluster_order).exists() else []
    if cats_json and not set(cats_json) & set(act.columns) and set(cats_json) & set(act.index):
        act = act.T
        act.columns = act.columns.map(str)
        act.index = act.index.map(str)

    cats = list(act.columns)
    cols = [c for c in cats_json if c in cats] + [c for c in cats if c not in cats_json] \
        if cats_json else sorted(cats, key=lambda s: int(s) if s.isdigit() else 1e9)

    # TF row order
    if args.tfs:
        want = list(dict.fromkeys(args.tfs))
        rows = [t for t in want if t in act.index]
        groups = {t: '' for t in rows}
        missing = [t for t in want if t not in act.index]
    else:
        rows, groups, missing = [], {}, []
        for grp, tfs in DEFAULT_TFS.items():
            for t in tfs:
                if t in act.index:
                    rows.append(t); groups[t] = grp
                else:
                    missing.append(t)
    if missing:
        print(f'[warn] TFs absent from the activity table (dropped): {missing}')
        print('       (if NKX2-1/ONECUT2/POU3F2 are missing, re-run decoupler_tf_activity.py '
              'with a prior that contains them, e.g. CollecTRI)')
    if not rows:
        raise SystemExit('no requested TFs present in the activity table')

    M = act.loc[rows, cols].astype(float)
    if args.scale == 'zscore':
        mu = M.mean(axis=1); sd = M.std(axis=1, ddof=0).replace(0, 1)
        D = M.sub(mu, axis=0).div(sd, axis=0)
        vmax = args.vlim; vmin = -args.vlim; cmap = 'RdBu_r'
        cbar = 'TF activity (z-score across clusters)'
        center = True
    else:
        D = M; vmax = float(np.nanmax(np.abs(D.values))); vmin = -vmax; cmap = 'RdBu_r'
        cbar = 'TF activity (decoupleR estimate)'; center = True

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt

    nrow, ncol = D.shape
    W = args.figwidth or (0.34 * ncol + 3.5)
    H = args.figheight or (0.26 * nrow + 1.8)
    fig, ax = plt.subplots(figsize=(W, H))
    im = ax.imshow(D.values, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(ncol)); ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(nrow)); ax.set_yticklabels(rows, fontsize=8)
    hl = {h.upper() for h in args.highlight}
    for t, r in zip(ax.get_yticklabels(), rows):
        if r.upper() in hl:
            t.set_fontweight('bold')
    colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if (args.cluster_colors and Path(args.cluster_colors).exists()) else {}
    if colors:
        for t, c in zip(ax.get_xticklabels(), cols):
            t.set_color(colors.get(c, 'black'))
    # program separators between TF groups (if grouped)
    if not args.tfs:
        seen, bnd = [], []
        for i, r in enumerate(rows):
            if groups[r] not in seen:
                seen.append(groups[r])
                if i > 0:
                    bnd.append(i)
        for b in bnd:
            ax.axhline(b - 0.5, color='0.5', lw=0.8)
        # group labels on the right
        tr = ax.get_yaxis_transform()
        starts = [0] + bnd + [nrow]
        for gi in range(len(seen)):
            mid = (starts[gi] + starts[gi + 1] - 1) / 2
            ax.text(1.005, mid, seen[gi], transform=tr, ha='left', va='center',
                    fontsize=7, rotation=90, color='0.3')

    ax.set_xlabel('cluster (adenocarcinoma → neuroendocrine)', fontsize=9)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.4, aspect=12, pad=0.10)
    cb.set_label(cbar, fontsize=8); cb.ax.tick_params(labelsize=7)
    ax.set_title('TF activity (decoupleR/CollecTRI): canonical NEPC drivers + MSX1',
                 fontsize=10)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}  ({nrow} TFs x {ncol} clusters)')


if __name__ == '__main__':
    main()
