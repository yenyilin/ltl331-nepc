#!/usr/bin/env python3
"""
compare_resolution_trajectory.py

Compare two (or more) Leiden resolutions for how a focus set of LTL331 original
clusters get grouped in the joint Harmony embedding. Default focus is the NEtD
trajectory: pre-EMT 4 and 6, the bridge 17, and early-NEPC 10 and 13.

Aggregate ARI tells you nothing about whether these five stay distinguishable.
This script answers, per resolution:
  1. Where does each focus cluster's LTL331 cells GO  (row-normalised crosstab).
  2. Do any two focus clusters collapse into the same joint cluster?
  3. For every joint cluster that absorbs focus cells, which ORIGINAL cluster
     dominates it across ALL LTL331 cells so cluster IDs are comparable across
     resolutions (a "dom=4" cluster at 0.5 is the analog of a "dom=4" cluster at
     0.55 even if the numeric IDs differ).

Outputs to OUTDIR (one set per resolution) plus a combined annotated landing
table for the focus clusters.

Usage:
    python compare_resolution_trajectory.py \\
        --h5ad integrated.h5ad \\
        --outdir resolution_compare \\
        --res-cols leiden_r0.5 leiden_r0.55 \\
        --group-col group --ltl331-label other \\
        --original-col original_cluster \\
        --focus 4 6 17 10 13
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def dominant_original_map(ltl_obs, original_col, res_col):
    """For each joint cluster, the dominant original cluster across ALL LTL331
    cells (not just the focus set). Used to label joint clusters comparably
    across resolutions."""
    ct = pd.crosstab(ltl_obs[res_col], ltl_obs[original_col])
    frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
    return pd.DataFrame({
        'dominant_original': frac.idxmax(axis=1),
        'dominant_frac': frac.max(axis=1).round(3),
        'n_ltl331_cells': ct.sum(axis=1).astype(int),
    })


def cohort_composition(obs, group_col, res_col, joint_clusters):
    """Cohort (Gao/Li/LTL331) make-up of given joint clusters, using ALL cells.

    Answers: of the joint clusters that LTL331 trajectory cells land in, what
    fraction of cells come from each cohort — i.e. do Gao/Li cells co-cluster
    there (cross-cohort validation)?
    """
    sub = obs[obs[res_col].isin(joint_clusters)]
    ct = pd.crosstab(sub[res_col], sub[group_col])
    ct = ct.reindex([c for c in joint_clusters if c in ct.index])
    ct['n_total'] = ct.sum(axis=1)
    cohort_cols = [c for c in ct.columns if c != 'n_total']
    pct = (ct[cohort_cols].div(ct['n_total'], axis=0) * 100).round(1)
    pct = pct.add_prefix('pct_')
    return pd.concat([ct, pct], axis=1)


def focus_landing(ltl_obs, original_col, res_col, focus, dom_map):
    """Row-normalised landing of each focus cluster, with destination joint
    clusters annotated by their dominant original cluster."""
    sub = ltl_obs[ltl_obs[original_col].isin(focus)]
    counts = pd.crosstab(sub[original_col], sub[res_col])
    counts = counts.reindex([f for f in focus if f in counts.index])
    counts = counts.loc[:, counts.sum(axis=0) > 0]

    rows = []
    for orig in counts.index:
        n_total = int(counts.loc[orig].sum())
        for joint in counts.columns:
            n = int(counts.loc[orig, joint])
            if n == 0:
                continue
            rows.append({
                'original_cluster': orig,
                'n_in_original': n_total,
                'joint_cluster': joint,
                'joint_dominant_original': dom_map.loc[joint, 'dominant_original'],
                'joint_dominant_frac': dom_map.loc[joint, 'dominant_frac'],
                'n_cells': n,
                'frac_of_original': round(n / n_total, 3),
            })
    landing = pd.DataFrame(rows).sort_values(
        ['original_cluster', 'frac_of_original'], ascending=[True, False]
    )
    return counts, landing


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__
    )
    p.add_argument('--h5ad', required=True)
    p.add_argument('--outdir', default='resolution_compare')
    p.add_argument('--res-cols', nargs='+', default=['leiden_r0.5', 'leiden_r0.55'],
                   help='obs columns with the clusterings to compare')
    p.add_argument('--group-col', default='group')
    p.add_argument('--ltl331-label', default='other',
                   help='value in --group-col identifying LTL331 cells')
    p.add_argument('--original-col', default='original_cluster')
    p.add_argument('--focus', nargs='+', default=['4', '6', '17', '10', '13'],
                   help='original cluster IDs to focus on (trajectory order)')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    obs[args.group_col] = obs[args.group_col].astype(str)
    obs[args.original_col] = obs[args.original_col].astype(str)
    for rc in args.res_cols:
        if rc not in obs.columns:
            raise ValueError(f'res column "{rc}" not in obs '
                             f'(have: {list(obs.columns)})')
        obs[rc] = obs[rc].astype(str)

    focus = [str(f) for f in args.focus]
    ltl = obs[obs[args.group_col] == args.ltl331_label]
    print(f'LTL331 cells (group=={args.ltl331_label}): {len(ltl)}')
    print(f'Focus original clusters: {focus}')
    for f in focus:
        n = int((ltl[args.original_col] == f).sum())
        print(f'  cluster {f}: n={n}')

    landings = {}
    for rc in args.res_cols:
        print('\n' + '=' * 70)
        print(f'RESOLUTION COLUMN: {rc}  '
              f'(n_joint_clusters={ltl[rc].nunique()})')
        print('=' * 70)

        dom = dominant_original_map(ltl, args.original_col, rc)
        dom.to_csv(outdir / f'dominant_original_{rc}.tsv', sep='\t')

        counts, landing = focus_landing(ltl, args.original_col, rc, focus, dom)
        counts.to_csv(outdir / f'focus_counts_{rc}.tsv', sep='\t')
        landing.to_csv(outdir / f'focus_landing_{rc}.tsv', sep='\t', index=False)
        landings[rc] = landing

        # Cohort composition of the destination joint clusters (ALL cells)
        dest_clusters = list(dict.fromkeys(landing['joint_cluster'].tolist()))
        comp = cohort_composition(obs, args.group_col, rc, dest_clusters)
        comp.insert(0, 'dom_original_ltl331',
                    dom.reindex(comp.index)['dominant_original'])
        comp.to_csv(outdir / f'cohort_composition_{rc}.tsv', sep='\t')

        # Pretty per-focus-cluster print
        for orig in focus:
            d = landing[landing['original_cluster'] == orig]
            if d.empty:
                continue
            n_total = int(d['n_in_original'].iloc[0])
            print(f'\noriginal cluster {orig} (n={n_total}) lands in:')
            for _, r in d.iterrows():
                tag = (f"joint {r['joint_cluster']} "
                       f"[dom={r['joint_dominant_original']} "
                       f"@{r['joint_dominant_frac']:.0%}]")
                print(f"    {r['frac_of_original']:>5.1%}  "
                      f"({r['n_cells']:>4d})  {tag}")

        # Cohort composition of destination joint clusters (cross-cohort view)
        print('\ncohort composition of destination joint clusters '
              '(all cells):')
        cohort_cols = [c for c in comp.columns if c.startswith('pct_')]
        for joint in comp.index:
            dtag = comp.loc[joint, 'dom_original_ltl331']
            n_tot = int(comp.loc[joint, 'n_total'])
            parts = '  '.join(
                f"{c.replace('pct_', '')} {comp.loc[joint, c]:.1f}%"
                for c in cohort_cols
            )
            print(f"    joint {joint} [dom={dtag}]  n={n_tot:>5d}   {parts}")

        # Collapse check: any joint cluster that is the top destination for >1
        # focus cluster?
        top_dest = (landing.sort_values('frac_of_original', ascending=False)
                            .groupby('original_cluster').first())
        dest_counts = top_dest['joint_cluster'].value_counts()
        merged = dest_counts[dest_counts > 1]
        if len(merged):
            print('\n[collapse warning] focus clusters sharing a top destination:')
            for joint, k in merged.items():
                who = top_dest[top_dest['joint_cluster'] == joint].index.tolist()
                dtag = dom.loc[joint, 'dominant_original']
                print(f'    joint {joint} (dom={dtag}) <- {who}')
        else:
            print('\n[ok] each focus cluster has a distinct top destination '
                  '(no collapse).')

    # Combined side-by-side landing for cluster 17 specifically (the bridge)
    bridge = focus[focus.index('17')] if '17' in focus else None
    if bridge:
        print('\n' + '=' * 70)
        print(f'SIDE-BY-SIDE: original cluster {bridge} landing across resolutions')
        print('=' * 70)
        for rc in args.res_cols:
            d = landings[rc]
            d17 = d[d['original_cluster'] == bridge]
            print(f'\n  {rc}:')
            for _, r in d17.iterrows():
                print(f"    {r['frac_of_original']:>5.1%}  ({r['n_cells']:>4d})  "
                      f"joint {r['joint_cluster']} "
                      f"[dom={r['joint_dominant_original']}]")

    print(f'\nSaved per-resolution TSVs to {outdir}/')
    print('Done.')


if __name__ == '__main__':
    main()
