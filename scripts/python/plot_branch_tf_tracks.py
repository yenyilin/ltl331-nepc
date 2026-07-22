#!/usr/bin/env python3
"""
plot_branch_tf_tracks.py — Fig 5C: TF activity along the two NE branches.

From the cluster-17 bridge, the trajectory commits at the cl13 NE-entry hub into
an ASCL1+ branch (cl13 -> cl10) and an ASCL1- branch (cl13 -> cl1 -> cl3 -> cl9
-> cl7) (see notes/ltl331_trajectory_model.md). This panel lays out decoupleR/
CollecTRI TF activity along each branch as two heatmaps that share the
bridge->entry trunk (cl17, cl13) and then diverge, so the reader sees which
regulators rise on which arm.

Input: tf_activity_per_cluster.tsv from decoupler_tf_activity.py (z-scored
pseudobulk recommended; rows = TF or cluster, auto-oriented). Color = per-TF
z-score across the UNION of branch clusters (one common scale for both panels so
the two arms are directly comparable). Canonical NE drivers + MSX1 bolded.

IMPORTANT (legend): the x-axis is CLUSTER STAGE along each branch, not continuous
per-cell pseudotime. It is the in-object summary of the branch; a continuous
per-cell pseudotime version requires the velocity object (cluster run).

Example:
  python plot_branch_tf_tracks.py \\
      --activity decoupler-results/all_methods_zscored/tf_activity_per_cluster.tsv \\
      --cluster-colors data/cluster_colors_18.json \\
      --outdir figures --stem fig5C_branch_tf_tracks
"""
import argparse, json
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

# branch definitions: ordered cluster stages (strings). Shared trunk = the
# clusters common to both, drawn at the left of each panel.
DEFAULT_BRANCHES = OrderedDict([
    ('ASCL1+ branch', ['17', '13', '10']),
    ('ASCL1− branch', ['17', '13', '1', '3', '9', '7']),
])
TRUNK = ['17', '13']  # bridge -> NE-entry, shared origin of both branches

# curated TFs grouped by program (PRAD -> EMT -> NE), plus the branch-
# discriminating neural-dev / engine TFs the run log surfaced (cl9 ATOH1/SIX1/
# POU4F3 neural core; cl3 E2F engine). Only those present are shown.
DEFAULT_TFS = OrderedDict([
    ('AR-axis / PRAD', ['AR', 'FOXA1', 'NKX3-1', 'GATA2']),
    ('EMT / bridge',   ['SNAI2', 'SNAI1', 'TWIST1', 'ZEB1']),
    ('NE drivers',     ['ASCL1', 'NEUROD1', 'INSM1', 'FOXA2', 'NKX2-1',
                        'POU3F2', 'SOX2', 'MSX1']),
    ('ASCL1− neural / engine', ['ATOH1', 'SIX1', 'POU4F3', 'E2F1']),
    ('NE repressor',   ['REST']),
])
HL = {'NKX2-1', 'POU3F2', 'ASCL1', 'NEUROD1', 'INSM1', 'FOXA2', 'MSX1'}


def load_activity(path, want_clusters):
    act = pd.read_csv(path, sep='\t', index_col=0)
    act.columns = act.columns.map(str)
    act.index = act.index.map(str)
    # orient so rows = TF, cols = cluster
    if set(want_clusters) & set(act.index) and not (set(want_clusters) & set(act.columns)):
        act = act.T
        act.columns = act.columns.map(str)
        act.index = act.index.map(str)
    return act


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--activity',
                    default='data/decoupler_tf/tf_activity_per_cluster.tsv')
    ap.add_argument('--cluster-colors', default='data/cluster_colors_18.json')
    ap.add_argument('--branches', default=None,
                    help='JSON file mapping branch name -> ordered list of cluster '
                         'ids; defaults to the cl17->cl13->{10}|{1,3,9,7} structure')
    ap.add_argument('--tfs', nargs='*', default=None,
                    help='explicit TF order (overrides the curated grouped default)')
    ap.add_argument('--vlim', type=float, default=2.0)
    ap.add_argument('--highlight', nargs='+', default=sorted(HL))
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig5C_branch_tf_tracks')
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

    if args.branches:
        branches = OrderedDict((k, [str(c) for c in v])
                               for k, v in json.load(open(args.branches)).items())
    else:
        branches = DEFAULT_BRANCHES

    all_clusters = list(dict.fromkeys(c for seq in branches.values() for c in seq))
    act = load_activity(args.activity, all_clusters)

    missing_cl = [c for c in all_clusters if c not in act.columns]
    if missing_cl:
        raise SystemExit(f'clusters not in activity table: {missing_cl}')

    # TF rows
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
    if not rows:
        raise SystemExit('no requested TFs present in the activity table')

    # z-score each TF across the UNION of branch clusters -> one common scale
    M = act.loc[rows, all_clusters].astype(float)
    mu = M.mean(axis=1); sd = M.std(axis=1, ddof=0).replace(0, 1)
    Z = M.sub(mu, axis=0).div(sd, axis=0)
    vmax = args.vlim; vmin = -args.vlim

    colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if (args.cluster_colors and Path(args.cluster_colors).exists()) else {}

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    bnames = list(branches.keys())
    nrow = len(rows)

    # ---- de-duplicated fork layout: the shared trunk (cl17->cl13) is drawn ONCE on the
    # left, then each arm forks off it (arm = branch minus the trunk). This removes the two
    # duplicated trunk columns that made the panel ~1/3 wider than it needs to be.
    trunk_seq = [c for c in TRUNK if c in all_clusters]
    arms = OrderedDict((b, [c for c in branches[b] if c not in TRUNK]) for b in bnames)
    # SHORT segment labels: the ASCL1+ arm is a single column, so a long title overlaps its
    # neighbours. Strip ' branch' -> 'ASCL1+' / 'ASCL1−'; trunk -> 'bridge'.
    seg_titles = ['bridge'] + [b.replace(' branch', '') for b in bnames]
    seg_seqs = [trunk_seq] + [arms[b] for b in bnames]
    widths = [max(len(s), 1) for s in seg_seqs]

    fig = plt.figure(figsize=figsize(0.40 * sum(widths) + 2.2, 0.30 * nrow + 1.6))
    gs = GridSpec(1, len(seg_seqs) + 1, width_ratios=widths + [0.35], wspace=0.12)

    im = None
    for si, (stitle, seq) in enumerate(zip(seg_titles, seg_seqs)):
        ax = fig.add_subplot(gs[0, si])
        D = Z.loc[rows, seq].values
        im = ax.imshow(D, aspect='auto', cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(seq)))
        ax.set_xticklabels(seq, fontsize=14)   # cluster IDs: enlarged for readability
        for t, c in zip(ax.get_xticklabels(), seq):
            t.set_color(colors.get(c, 'black'))
            if c in TRUNK:
                t.set_style('italic')
        ax.set_title(stitle, fontsize=8)
        if si == 0:
            ax.set_yticks(range(nrow)); ax.set_yticklabels(rows, fontsize=8)
            hl = {h.upper() for h in args.highlight}
            for t, r in zip(ax.get_yticklabels(), rows):
                if r.upper() in hl:
                    t.set_fontweight('bold')
            # program separators (only on leftmost / trunk panel)
            if not args.tfs:
                seen, bnd = [], []
                for i, r in enumerate(rows):
                    if groups[r] not in seen:
                        seen.append(groups[r])
                        if i > 0:
                            bnd.append(i)
                for b in bnd:
                    ax.axhline(b - 0.5, color='0.5', lw=0.8)
        else:
            ax.set_yticks([])
        ax.tick_params(length=0)

    cax = fig.add_subplot(gs[0, -1])
    cb = fig.colorbar(im, cax=cax, fraction=1.0)
    cb.set_label('TF activity (z)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    # squeeze the colorbar axis to a slim bar
    pos = cax.get_position()
    cax.set_position([pos.x0, pos.y0 + pos.height * 0.25,
                      pos.width * 0.35, pos.height * 0.5])

    fig.supxlabel('cluster stage', fontsize=8)
    # short title so bbox_inches='tight' can't hold the panel wide (the old one-liner did)
    fig.suptitle('TF activity along the NE branches', fontsize=9, y=1.01)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}  '
          f'({nrow} TFs; branches: ' +
          '; '.join(f'{b}={">".join(branches[b])}' for b in bnames) + ')')


if __name__ == '__main__':
    main()
