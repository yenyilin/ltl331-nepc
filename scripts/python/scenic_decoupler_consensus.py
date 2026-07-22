#!/usr/bin/env python3
"""
scenic_decoupler_consensus.py — cross-method consensus for cluster-specific
TF regulators. Cross-validates SCENIC (co-expression motif enrichment) against
decoupleR (prior-knowledge TF activity) to identify high-confidence cluster
regulators.

For each (cluster, TF) pair, a status is assigned:
  CONCORDANT      SCENIC rank ≤ --scenic-topk AND decoupleR |z| ≥ --decoupler-zthresh
                  (and decoupleR |z| rank within cluster ≤ --decoupler-topk),
                  WITH consistent direction (SCENIC `_(+)` matches z>0; `_(-)`
                  matches z<0). Bulletproof "active in both methods".
  DIVERGENT       Both methods flag the TF as cluster-specific BUT directions
                  disagree. Worth a closer look — possible prior-network sign
                  error, or TF doing something unusual in this dataset.
  SCENIC-only     SCENIC top-K but no decoupleR signal. Common for homeobox /
                  developmental TFs with sparse CollecTRI regulons.
  decoupleR-only  Active in z-scored ULM but not in SCENIC top-K. Indicates a
                  prior-knowledge hit not picked up by co-expression motifs.

SCENIC regulons of the form `TF_(+)` and `TF_(-)` for the same TF are
collapsed by taking the most cluster-specific form (lowest rank).

Inputs
  --scenic-rank   data/scenic/all_all_rss_rank.tsv (regulons × clusters,
                  rank values; regulon names suffixed _(+) / _(-))
  --decoupler     decoupler-results/ulm_zscored/tf_activity_per_cluster.tsv
                  (clusters × TFs, signed z-scores from decoupleR ULM on
                  z-scored pseudobulk)

Outputs (in --outdir)
  consensus_per_cluster.tsv       long: cluster, tf, status, scenic_rank,
                                  regulon_dir, decoupler_z, decoupler_rank_abs,
                                  regulon_full
  consensus_summary.tsv           per-cluster counts of each status
  high_confidence_regulators.tsv  CONCORDANT-only rows, the per-cluster list
                                  for future inspection
  consensus_canonical_tfs.tsv     wide: canonical NEPC TFs × clusters with
                                  (scenic_rank, decoupler_z, status) columns

Example
  python3 scripts/scenic_decoupler_consensus.py \\
      --scenic-rank data/scenic/all_all_rss_rank.tsv \\
      --decoupler decoupler-results/ulm_zscored/tf_activity_per_cluster.tsv \\
      --scenic-topk 10 --decoupler-topk 20 --decoupler-zthresh 1.5 \\
      --outdir consensus_tf
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


CANONICAL_NEPC_TFS = [
    'ASCL1', 'NEUROD1', 'POU3F2', 'POU2F3', 'INSM1', 'FOXA2', 'NKX2-1',
    'ONECUT2', 'REST', 'SOX2', 'MSX1', 'YAP1',
    'AR', 'FOXA1', 'NKX3-1', 'TP63', 'TP53', 'MYC',
]


def parse_regulon(name):
    """Parse a SCENIC regulon name like 'ASCL1_(+)' -> ('ASCL1', '+')."""
    m = re.match(r'^(.+?)_\(([+-])\)$', name)
    if not m:
        return name, ''
    return m.group(1), m.group(2)


def numkey(c):
    """Sort key: numeric-looking strings first, sorted numerically."""
    s = str(c)
    return (0, int(s)) if s.isdigit() else (1, s)


def load_scenic_rank(path):
    """Long DataFrame: cluster, tf, regulon_dir, scenic_rank, regulon_full.
    Both `_(+)` and `_(-)` regulons for the same TF are kept as separate rows;
    collapsed downstream by taking the most cluster-specific (lowest rank).
    """
    df = pd.read_csv(path, sep='\t', index_col=0)
    df.columns = df.columns.astype(str)
    rows = []
    for reg in df.index:
        tf, sign = parse_regulon(reg)
        for cl in df.columns:
            rows.append({'cluster': cl, 'tf': tf, 'regulon_dir': sign,
                         'scenic_rank': float(df.loc[reg, cl]),
                         'regulon_full': reg})
    return pd.DataFrame(rows)


def load_decoupler(path):
    """Long DataFrame: cluster, tf, decoupler_z, decoupler_rank_abs.
    decoupler_rank_abs is rank within cluster by |z| (1 = most cluster-specific).
    """
    df = pd.read_csv(path, sep='\t', index_col=0)
    df.index = df.index.astype(str)
    rows = []
    for cl in df.index:
        z = df.loc[cl]
        ranks = z.abs().rank(ascending=False, method='min').astype(int)
        for tf in df.columns:
            rows.append({'cluster': cl, 'tf': tf,
                         'decoupler_z': float(z[tf]),
                         'decoupler_rank_abs': int(ranks[tf])})
    return pd.DataFrame(rows)


def collapse_scenic_per_tf(scenic_long):
    """Keep the most cluster-specific regulon (lowest rank) per (cluster, tf)."""
    idx = scenic_long.groupby(['cluster', 'tf'])['scenic_rank'].idxmin()
    return scenic_long.loc[idx].reset_index(drop=True)


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--scenic-rank', required=True,
                   help='SCENIC RSS rank TSV (regulons × clusters)')
    p.add_argument('--decoupler', required=True,
                   help='decoupleR z-scored activity TSV (clusters × TFs)')
    p.add_argument('--scenic-topk', type=int, default=10,
                   help='SCENIC cluster-specific cutoff (rank ≤ this; default 10)')
    p.add_argument('--decoupler-topk', type=int, default=20,
                   help='decoupleR cluster-specific cutoff (|z| rank ≤ this; '
                        'default 20)')
    p.add_argument('--decoupler-zthresh', type=float, default=1.5,
                   help='also require |decoupler z| ≥ this for CONCORDANT '
                        '(default 1.5)')
    p.add_argument('--canonical-tfs', nargs='+', default=CANONICAL_NEPC_TFS)
    p.add_argument('--outdir', default='consensus_tf')
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # 1. load
    scenic = load_scenic_rank(args.scenic_rank)
    deco = load_decoupler(args.decoupler)
    print(f'SCENIC : {scenic["tf"].nunique()} TFs, '
          f'{scenic["cluster"].nunique()} clusters, '
          f'{scenic.shape[0]:,} (cluster, regulon) rows')
    print(f'decoupleR: {deco["tf"].nunique()} TFs, '
          f'{deco["cluster"].nunique()} clusters')

    # 2. collapse SCENIC per (cluster, TF) to the most-specific regulon
    scenic_collapsed = collapse_scenic_per_tf(scenic)

    # 3. outer-merge so TFs in only one method are kept
    merged = scenic_collapsed.merge(deco, on=['cluster', 'tf'], how='outer')

    # 4. tag each row
    sk = args.scenic_topk
    dk = args.decoupler_topk
    dz = args.decoupler_zthresh

    def tag(row):
        in_scenic = pd.notna(row['scenic_rank']) and row['scenic_rank'] <= sk
        in_deco = (pd.notna(row['decoupler_rank_abs'])
                   and row['decoupler_rank_abs'] <= dk
                   and pd.notna(row['decoupler_z'])
                   and abs(row['decoupler_z']) >= dz)
        if in_scenic and in_deco:
            scenic_pos = (row['regulon_dir'] == '+')
            deco_pos = row['decoupler_z'] > 0
            return 'CONCORDANT' if scenic_pos == deco_pos else 'DIVERGENT'
        if in_scenic:
            return 'SCENIC-only'
        if in_deco:
            return 'decoupleR-only'
        return 'neither'

    merged['status'] = merged.apply(tag, axis=1)
    merged = merged.sort_values(
        ['cluster', 'status', 'scenic_rank'], na_position='last')

    # 5. main long output (exclude "neither")
    main_out = merged[merged['status'] != 'neither'].copy()
    cols = ['cluster', 'tf', 'status', 'scenic_rank', 'regulon_dir',
            'decoupler_z', 'decoupler_rank_abs', 'regulon_full']
    cols = [c for c in cols if c in main_out.columns]
    main_out[cols].to_csv(outdir / 'consensus_per_cluster.tsv',
                          sep='\t', index=False)
    print(f'\nwrote {outdir}/consensus_per_cluster.tsv '
          f'({len(main_out):,} non-neither rows)')

    # 6. per-cluster status counts
    summary = (main_out.groupby(['cluster', 'status']).size()
               .unstack(fill_value=0)
               .reindex(columns=['CONCORDANT', 'DIVERGENT', 'SCENIC-only',
                                 'decoupleR-only'], fill_value=0))
    summary = summary.reindex(sorted(summary.index, key=numkey))
    summary.to_csv(outdir / 'consensus_summary.tsv', sep='\t')
    print('\nper-cluster status counts:')
    print(summary.to_string())

    # 7. high-confidence regulators (CONCORDANT only)
    hc = main_out[main_out['status'] == 'CONCORDANT'].copy()
    hc = hc.sort_values(['cluster', 'scenic_rank'],
                        key=lambda s: (s.map(numkey) if s.name == 'cluster'
                                       else s))
    hc[cols].to_csv(outdir / 'high_confidence_regulators.tsv',
                    sep='\t', index=False)
    print(f'\nwrote {outdir}/high_confidence_regulators.tsv '
          f'({len(hc):,} CONCORDANT TF-cluster pairs)')

    # 8. canonical-NEPC-TF wide cross-method table
    canset = set(args.canonical_tfs)
    cano = merged[merged['tf'].isin(canset)].copy()
    clusters_sorted = sorted(merged['cluster'].dropna().unique(), key=numkey)
    wide_rows = []
    for tf in args.canonical_tfs:
        row = {'TF': tf}
        for cl in clusters_sorted:
            sub = cano[(cano['tf'] == tf) & (cano['cluster'] == cl)]
            if sub.empty:
                row[f'cl{cl}_scenic_rank'] = np.nan
                row[f'cl{cl}_decoupler_z'] = np.nan
                row[f'cl{cl}_status'] = ''
            else:
                r = sub.iloc[0]
                row[f'cl{cl}_scenic_rank'] = (
                    int(r['scenic_rank']) if pd.notna(r['scenic_rank']) else np.nan)
                row[f'cl{cl}_decoupler_z'] = (
                    round(float(r['decoupler_z']), 2)
                    if pd.notna(r['decoupler_z']) else np.nan)
                row[f'cl{cl}_status'] = r['status'] if r['status'] != 'neither' else ''
        wide_rows.append(row)
    wide = pd.DataFrame(wide_rows)
    wide.to_csv(outdir / 'consensus_canonical_tfs.tsv', sep='\t', index=False)
    print(f'wrote {outdir}/consensus_canonical_tfs.tsv')

    # 9. printed summary of canonical NEPC TF concordance
    print('\n=== canonical NEPC TF concordance ===')
    cano_hc = cano[cano['status'] == 'CONCORDANT']
    cano_div = cano[cano['status'] == 'DIVERGENT']
    cano_scenic = cano[cano['status'] == 'SCENIC-only']
    cano_deco = cano[cano['status'] == 'decoupleR-only']

    if not cano_hc.empty:
        print('\nCONCORDANT (both methods, same direction):')
        print(cano_hc[['tf', 'cluster', 'scenic_rank', 'regulon_dir',
                       'decoupler_z']].to_string(index=False))
    if not cano_div.empty:
        print('\nDIVERGENT (both methods, opposite direction):')
        print(cano_div[['tf', 'cluster', 'scenic_rank', 'regulon_dir',
                        'decoupler_z']].to_string(index=False))
    if not cano_scenic.empty:
        print('\nSCENIC-only (no decoupleR corroboration):')
        print(cano_scenic[['tf', 'cluster', 'scenic_rank',
                           'regulon_dir']].to_string(index=False))
    if not cano_deco.empty:
        print('\ndecoupleR-only (no SCENIC corroboration):')
        print(cano_deco[['tf', 'cluster', 'decoupler_z',
                         'decoupler_rank_abs']].to_string(index=False))

    flagged = (cano[cano['status'] != 'neither']['tf'].unique())
    never_hit = [tf for tf in args.canonical_tfs if tf not in flagged]
    if never_hit:
        print(f'\ncanonical NEPC TFs NOT flagged in ANY cluster by either method: '
              f'{never_hit}')


if __name__ == '__main__':
    main()
