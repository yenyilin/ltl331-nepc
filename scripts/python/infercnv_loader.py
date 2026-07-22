#!/usr/bin/env python3
"""
infercnv_loader.py — parse an inferCNV HMM output directory and emit a tidy
long-format CNV calls TSV matching the schema of nexus_cnv_loader.py, so the
existing plot_nexus_cnv_heatmap.py / plot_nexus_cnv_curves.py scripts work on
inferCNV calls without modification.

Inputs (under --infercnv-dir, file naming as inferCNV ships them):
  HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat
      per-subcluster CNV segments (chr, start, end, HMMi6 state 1..6)
  17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings
      cell -> subcluster mapping (used for cell-count weighting in
      per-timepoint aggregation)

HMMi6 state encoding -> Nexus-style integer codes:
  1 = complete loss     -> -2 (big_loss)
  2 = 1-copy loss       -> -1 (loss)
  3 = neutral           ->  0  (filtered out of segment output)
  4 = 1-copy gain       -> +1 (gain)
  5 = 2-copy gain       -> +2 (high_gain)
  6 = multi-copy        -> +2 (high_gain, capped)

cell_group_name in inferCNV looks like "preCx1.preCx1_s2"; the group prefix is
the metadata column passed in (here: timepoint), the suffix `s\d+` is the
inferCNV subcluster id.  Timepoint strings are normalised (preCx1 -> preCx,
16wk1 -> 16wk) so they line up with the WES sample names emitted by
nexus_cnv_loader.py.

Aggregation modes (--by):
  timepoint (default)
      Pseudo-bulk per timepoint.  For each bin (default 1 Mbp) the subcluster
      CNV codes are averaged weighted by cell count, then rounded to
      {-2,-1,0,+1,+2}, then run-length-encoded back to segments.  Subclusters
      with no entry for a bin are treated as neutral (0).
  subcluster
      Raw per-subcluster calls passed through.  Each subcluster becomes one
      "sample" in the output TSV (named "<timepoint>.s<n>").

Outputs:
  <prefix>_calls.tsv         long-format (sample, Chromosome, Start, End, Value, value_label)
  <prefix>_descriptors.tsv   one row per sample (timepoint or subcluster) with
                              parsed metadata (group, n_cells, ...)

Example:
  python infercnv_loader.py \\
      --infercnv-dir data/infercnv/samples.groups \\
      --out-prefix data/infercnv \\
      --by timepoint --bin-size 1000000
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from plot_nexus_cnv_heatmap import (
    HG19_SIZES, CHROM_ORDER, build_bins, paint_sample,
)


REGIONS_FNAME = ('HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters'
                 '.Pnorm_0.5.pred_cnv_regions.dat')
CELLS_FNAME = '17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings'

HMM_TO_NEXUS = {1: -2, 2: -1, 3: 0, 4: 1, 5: 2, 6: 2}
LABEL = {-2: 'big_loss', -1: 'loss', 0: 'neutral', 1: 'gain', 2: 'high_gain'}


# --------------------------------------------------------------------------- #
# parsing helpers
# --------------------------------------------------------------------------- #
def _parse_cell_group(s):
    """'preCx1.preCx1_s2' -> ('preCx1', 's2')."""
    if not isinstance(s, str) or '.' not in s:
        return s, None
    group, sub = s.split('.', 1)
    m = re.search(r'_(s\d+)$', sub)
    return group, m.group(1) if m else sub


def _normalize_timepoint(tp):
    """Harmonise inferCNV group names with the WES timepoint convention used
    by nexus_cnv_loader.py.

      preCx, preCX1, preCx2 -> preCx
      week08, week8R        -> 8wk           (also strips leading zeros + R suffix)
      week20R, week22R      -> 20wk, 22wk    (R = regression/relapse cohort;
                                              same biological timepoint)
      16wk1                  -> 16wk
      d17                    -> d17           (no WES counterpart; left as-is)
    """
    if not isinstance(tp, str):
        return tp
    t = tp
    t = re.sub(r'preCX?\d*', 'preCx', t, flags=re.IGNORECASE)
    m = re.match(r'^week0*(\d+)R?$', t, flags=re.IGNORECASE)
    if m:
        t = f'{m.group(1)}wk'
    t = re.sub(r'(\d+wk)\d+$', r'\1', t)
    return t


def _runs(vec):
    """RLE: yield (start_idx, end_idx_inclusive, value) for consecutive runs."""
    n = len(vec)
    if n == 0:
        return
    start = 0
    for i in range(1, n):
        if vec[i] != vec[start]:
            yield (start, i - 1, int(vec[start]))
            start = i
    yield (start, n - 1, int(vec[start]))


# --------------------------------------------------------------------------- #
# loaders
# --------------------------------------------------------------------------- #
def load_inputs(infercnv_dir):
    d = Path(infercnv_dir)
    regions = pd.read_csv(d / REGIONS_FNAME, sep='\t')
    cells = pd.read_csv(d / CELLS_FNAME, sep='\t')
    return regions, cells


def recode(regions):
    out = regions.copy()
    out['Value'] = out['state'].map(HMM_TO_NEXUS).astype(int)
    out['value_label'] = out['Value'].map(LABEL)
    parsed = out['cell_group_name'].apply(_parse_cell_group)
    out['group_raw'] = parsed.apply(lambda x: x[0])
    out['subcluster'] = parsed.apply(lambda x: x[1])
    out['timepoint'] = out['group_raw'].apply(_normalize_timepoint)
    return out


def cell_counts(cells):
    parsed = cells['cell_group_name'].apply(_parse_cell_group)
    cells = cells.copy()
    cells['group_raw'] = parsed.apply(lambda x: x[0])
    cells['subcluster'] = parsed.apply(lambda x: x[1])
    cells['timepoint'] = cells['group_raw'].apply(_normalize_timepoint)
    return (cells.groupby(['timepoint', 'subcluster']).size()
                  .rename('n_cells').reset_index())


# --------------------------------------------------------------------------- #
# aggregation modes
# --------------------------------------------------------------------------- #
def per_subcluster_calls(regions):
    """One row per (subcluster, region); raw recoded calls, neutrals dropped."""
    df = regions[regions['Value'] != 0].copy()
    df['sample'] = df['timepoint'].astype(str) + '.' + df['subcluster'].fillna('s?')
    out = df[['sample', 'chr', 'start', 'end', 'Value', 'value_label']].rename(
        columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})
    return out


def aggregate_per_timepoint(regions, counts, bin_size):
    """Per-timepoint weighted-mean of subcluster CNV codes, RLE'd to segments.

    For each timepoint:
      - paint each subcluster's regions into a per-bin int8 vector
      - weighted average across subclusters (weights = cell counts)
      - round to {-2,-1,0,+1,+2}
      - RLE within each chromosome -> emit non-neutral segments
    """
    bins_df = build_bins(bin_size)
    n_total = len(bins_df)
    chrom_offsets, chrom_nbins = {}, {}
    for chrom, group in bins_df.groupby('chrom', sort=False):
        chrom_offsets[chrom] = int(group['global_idx'].min())
        chrom_nbins[chrom] = len(group)

    rows = []
    for tp, tp_regions in regions.groupby('timepoint'):
        tp_subs_in_regions = set(tp_regions['subcluster'].dropna().unique())
        tp_subs_in_cells = set(counts[counts['timepoint'] == tp]['subcluster'])
        all_subs = sorted(tp_subs_in_regions | tp_subs_in_cells)
        if not all_subs:
            continue

        sub_matrix = np.zeros((len(all_subs), n_total), dtype=np.float32)
        weights = np.zeros(len(all_subs), dtype=np.float32)
        for i, sub in enumerate(all_subs):
            sub_regs = tp_regions[tp_regions['subcluster'] == sub][
                ['chr', 'start', 'end', 'Value']
            ].rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})
            if not sub_regs.empty:
                painted = paint_sample(sub_regs, bin_size,
                                       chrom_offsets, chrom_nbins, n_total)
                sub_matrix[i] = painted
            cc = counts[(counts['timepoint'] == tp) &
                        (counts['subcluster'] == sub)]['n_cells']
            weights[i] = float(cc.iloc[0]) if len(cc) else 0.0

        if weights.sum() == 0:
            continue
        weighted = (sub_matrix * weights[:, None]).sum(axis=0) / weights.sum()
        rounded = np.round(weighted).astype(int).clip(-2, 2)

        # RLE per chromosome -> segments
        for chrom in CHROM_ORDER:
            cmask = (bins_df['chrom'] == chrom).values
            if not cmask.any():
                continue
            cidx_offset = int(np.where(cmask)[0].min())
            cvec = rounded[cmask]
            for s, e, v in _runs(cvec):
                if v == 0:
                    continue
                rows.append({
                    'sample': tp, 'Chromosome': chrom,
                    'Start': int(bins_df.iloc[cidx_offset + s]['start']),
                    'End': int(bins_df.iloc[cidx_offset + e]['end']),
                    'Value': int(v), 'value_label': LABEL[int(v)],
                })

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# descriptors
# --------------------------------------------------------------------------- #
def make_descriptors(calls, counts, by):
    samples = sorted(calls['sample'].unique())
    rows = []
    if by == 'timepoint':
        per_tp = counts.groupby('timepoint')['n_cells'].sum()
        for tp in samples:
            rows.append({'sample': tp, 'timepoint': tp,
                         'subcluster': None,
                         'n_cells': int(per_tp.get(tp, 0))})
    else:  # subcluster
        for s in samples:
            tp, sub = (s.split('.', 1) + [None])[:2]
            n = counts[(counts['timepoint'] == tp) &
                       (counts['subcluster'] == sub)]['n_cells']
            rows.append({'sample': s, 'timepoint': tp,
                         'subcluster': sub,
                         'n_cells': int(n.iloc[0]) if len(n) else None})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--infercnv-dir', required=True,
                   help='inferCNV HMM output directory')
    p.add_argument('--out-prefix', required=True,
                   help='write <prefix>_calls.tsv / _descriptors.tsv')
    p.add_argument('--by', choices=['timepoint', 'subcluster'], default='timepoint',
                   help='aggregation: per-timepoint pseudo-bulk (default) or '
                        'raw per-subcluster')
    p.add_argument('--bin-size', type=int, default=1_000_000,
                   help='bin width for per-timepoint aggregation (default 1 Mbp; '
                        'ignored when --by subcluster)')
    args = p.parse_args()

    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    regions, cells = load_inputs(args.infercnv_dir)
    print(f'Loaded {len(regions)} subcluster CNV regions, {len(cells)} cells')
    regions = recode(regions)
    counts = cell_counts(cells)

    print('\nHMMi6 state -> Nexus code distribution:')
    print(regions['value_label'].value_counts().to_string())
    print(f'\nGroups (raw -> normalised timepoint):')
    print(regions[['group_raw', 'timepoint']].drop_duplicates().to_string(index=False))

    if args.by == 'subcluster':
        calls = per_subcluster_calls(regions)
    else:
        calls = aggregate_per_timepoint(regions, counts, args.bin_size)
    descr = make_descriptors(calls, counts, args.by)

    calls.to_csv(f'{args.out_prefix}_calls.tsv', sep='\t', index=False)
    descr.to_csv(f'{args.out_prefix}_descriptors.tsv', sep='\t', index=False)
    print(f'\nWrote {args.out_prefix}_calls.tsv '
          f'({len(calls)} segments, {calls["sample"].nunique()} samples)')
    print(f'Wrote {args.out_prefix}_descriptors.tsv')

    if not calls.empty:
        print('\nPer-sample segment count:')
        print(calls.groupby('sample').size().to_string())


if __name__ == '__main__':
    main()
