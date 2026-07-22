#!/usr/bin/env python3
"""
nexus_cnv_loader.py — walk a Nexus BioDiscovery export and assemble per-sample
CNV calls + zygosity (LOH) + sample metadata into tidy long-format DataFrames
ready for custom plotting (ideograms, arm-level frequency, sample-x-genome
heatmaps, LOH overlays).

Each `Samples/<name>/` folder is expected to contain:
  calls.txt              segment CNV calls    (Chromosome, Start, End, Value)
  zygosity.txt           LOH / allelic imbalance segments (same columns)
  descriptor.set         per-sample metadata (MAPD, Quality, Gender, BAM paths)
  callthresholds.txt     per-chromosome thresholds used for the calls (not loaded
                         by default; pass --keep-thresholds to include)
  *.ivg                  proprietary Nexus binary (per-probe log2R + BAF) — not
                         readable without Nexus and intentionally ignored here.

CNV `Value` codes (Nexus SNP-FASST2 default):
   -2 Big Loss     -1 Loss     +1 Gain     +2 High Gain

Zygosity `Value` codes:
   0 Homozygous Low     1 Homozygous High     3 Allelic Imbalance (LOH proxy)

Reference build per the manuscript: hg38 / NCBI Build 38 (see each sample's
settings.xml for confirmation).

Usable two ways:

  (a) CLI — write three TSVs for downstream plotting scripts:
      python nexus_cnv_loader.py \\
          --nexus-dir nexus/927-331.WGS \\
          --out-prefix data/nexus \\
          --filter-prefix LTL331

  (b) Import — load DataFrames directly:
      from nexus_cnv_loader import walk_samples
      calls, zyg, descr = walk_samples(
          'nexus/927-331.WGS', filter_prefix='LTL331')
"""
import argparse
import re
from pathlib import Path

import pandas as pd


CALL_LABEL = {-2: 'big_loss', -1: 'loss', 1: 'gain', 2: 'high_gain'}
ZYG_LABEL = {0: 'hom_low', 1: 'hom_high', 3: 'imbalance'}


# --------------------------------------------------------------------------- #
# parsing helpers
# --------------------------------------------------------------------------- #
def _read_segments(path, value_label_map):
    df = pd.read_csv(path, sep='\t')
    df['value_label'] = df['Value'].map(value_label_map)
    return df


def _read_descriptor(path):
    """descriptor.set is `Data Type:<tab>...` then a header row then one data
    row. Returns a dict; empty dict if the file is malformed/empty.
    """
    try:
        df = pd.read_csv(path, sep='\t', skiprows=1)
    except Exception:
        return {}
    if df.empty:
        return {}
    return df.iloc[0].to_dict()


def _parse_sample_name(folder_name):
    """Split a Nexus folder name like 'LTL331_preCX1 vs <matched normal>' into fields.

    Best-effort: non-LTL331 names yield NaN timepoint.  Recognised LTL331
    timepoint tokens follow the manuscript convention (preCx, 8wk, 12wk, 16wk,
    20wk, 22wk).  Trailing replicate digits (preCX1 -> preCx, 16wk1 -> 16wk)
    are stripped so the column lines up cleanly with timepoint groupings.
    """
    m = re.match(r'(.+?)\s+vs\s+(.+)', folder_name)
    if m:
        tumor, normal = m.group(1).strip(), m.group(2).strip()
    else:
        tumor, normal = folder_name, None

    series, timepoint = None, None
    if 'LTL331' in tumor:
        series = 'LTL331'
        tp = re.search(r'(preCX?\d*|\d+wk\d*)', tumor, flags=re.IGNORECASE)
        if tp:
            t = tp.group(1)
            t = re.sub(r'preCX?\d*', 'preCx', t, flags=re.IGNORECASE)
            t = re.sub(r'(\d+wk)\d+$', r'\1', t)
            timepoint = t
    elif tumor.startswith('927'):
        series = '927'
        # "927 WGS" is the parent prostate tumor (pre-castration baseline that
        # gave rise to the long-term line from which the LTL331 PDX series derives).
        # Tag it 'parent' so sort_samples can place it before preCx.
        if tumor.lower().startswith('927 wgs'):
            timepoint = 'parent'
    elif tumor.startswith('331-'):
        series = '331'

    # public/clean sample name = tumor only (drop the " vs <reference>" suffix so
    # the normalization-reference name is never written to hosted files)
    sample = '927' if tumor.strip().lower() == '927 wgs' else tumor.strip()
    return {
        'sample': sample,
        'sample_raw': folder_name,
        'tumor': tumor, 'normal': normal,
        'series': series, 'timepoint': timepoint,
    }


# --------------------------------------------------------------------------- #
# main walker
# --------------------------------------------------------------------------- #
def walk_samples(nexus_dir, filter_prefix=None, keep_thresholds=False):
    """Walk `<nexus_dir>/Samples/*/` and return (calls, zyg, descr).

    Parameters
    ----------
    nexus_dir : str | Path
        Nexus export root (the folder that contains a `Samples/` subdir).
    filter_prefix : str | list[str] | None
        If given, only sample folders whose name starts with any of these
        prefixes are loaded (e.g. ['LTL331', '927 WGS']).
    keep_thresholds : bool
        Also concatenate `callthresholds.txt` into a 4th DataFrame
        (returned when True; otherwise omitted).

    Returns
    -------
    (calls, zyg, descr)  or  (calls, zyg, descr, thresholds) if
    keep_thresholds=True.
    """
    samples_dir = Path(nexus_dir) / 'Samples'
    if not samples_dir.exists():
        raise SystemExit(f'no Samples/ found under {nexus_dir}')

    if isinstance(filter_prefix, str):
        filter_prefix = [filter_prefix]

    calls_rows, zyg_rows, descr_rows, thr_rows = [], [], [], []
    for sample_dir in sorted(samples_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        name = sample_dir.name
        if filter_prefix and not any(name.startswith(p) for p in filter_prefix):
            continue

        parsed = _parse_sample_name(name)
        clean = parsed['sample']   # public name, no " vs <reference>" suffix

        calls_path = sample_dir / 'calls.txt'
        if calls_path.exists():
            df = _read_segments(calls_path, CALL_LABEL)
            df['sample'] = clean
            calls_rows.append(df)

        zyg_path = sample_dir / 'zygosity.txt'
        if zyg_path.exists():
            df = _read_segments(zyg_path, ZYG_LABEL)
            df['sample'] = clean
            zyg_rows.append(df)

        descr_path = sample_dir / 'descriptor.set'
        if descr_path.exists():
            descr = _read_descriptor(descr_path)
            descr.update(parsed)   # sets clean 'sample' + series/timepoint
            descr_rows.append(descr)

        if keep_thresholds:
            thr_path = sample_dir / 'callthresholds.txt'
            if thr_path.exists():
                df = pd.read_csv(thr_path, sep='\t')
                df['sample'] = clean
                thr_rows.append(df)

    calls = pd.concat(calls_rows, ignore_index=True) if calls_rows else pd.DataFrame()
    zyg = pd.concat(zyg_rows, ignore_index=True) if zyg_rows else pd.DataFrame()
    descr = pd.DataFrame(descr_rows)
    # public-safe descriptor: drop the columns that reveal the normalization
    # reference or internal file paths (normalization-reference names, BAM paths, "X vs Y" ids)
    _sensitive = ['sample_raw', 'tumor', 'normal', 'Sample Name', 'File', 'Matched File']
    descr = descr.drop(columns=[c for c in _sensitive if c in descr.columns],
                       errors='ignore')

    cols_order = ['sample', 'Chromosome', 'Start', 'End', 'Value', 'value_label']
    if not calls.empty:
        calls = calls[cols_order]
    if not zyg.empty:
        zyg = zyg[cols_order]

    if keep_thresholds:
        thr = pd.concat(thr_rows, ignore_index=True) if thr_rows else pd.DataFrame()
        return calls, zyg, descr, thr
    return calls, zyg, descr


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--nexus-dir', required=True,
                   help='Nexus export root (folder containing Samples/)')
    p.add_argument('--out-prefix', default=None,
                   help='write <prefix>_calls.tsv / _zygosity.tsv / '
                        '_descriptors.tsv (omit to print summary only)')
    p.add_argument('--filter-prefix', nargs='+', default=None,
                   help='only sample folders starting with ANY of these '
                        'prefixes (e.g. LTL331 "927 WGS")')
    p.add_argument('--keep-thresholds', action='store_true',
                   help='also concatenate per-chromosome callthresholds.txt')
    args = p.parse_args()

    result = walk_samples(args.nexus_dir,
                          filter_prefix=args.filter_prefix,
                          keep_thresholds=args.keep_thresholds)
    if args.keep_thresholds:
        calls, zyg, descr, thr = result
    else:
        calls, zyg, descr = result
        thr = None

    print('Loaded:')
    print(f'  calls:       {len(calls):>6d} segments across '
          f'{calls["sample"].nunique() if not calls.empty else 0} samples')
    print(f'  zygosity:    {len(zyg):>6d} segments across '
          f'{zyg["sample"].nunique() if not zyg.empty else 0} samples')
    print(f'  descriptors: {len(descr):>6d} samples')
    if thr is not None:
        print(f'  thresholds:  {len(thr):>6d} chrom-rows')

    if not calls.empty:
        print('\nCall-value distribution (across all loaded samples):')
        print(calls['value_label'].value_counts().to_string())

    if not descr.empty:
        keep = [c for c in ('sample', 'series', 'timepoint', 'tumor', 'normal')
                if c in descr.columns]
        print('\nParsed sample identifiers:')
        print(descr[keep].to_string(index=False))

    if args.out_prefix:
        out = Path(args.out_prefix)
        out.parent.mkdir(parents=True, exist_ok=True)
        calls.to_csv(f'{args.out_prefix}_calls.tsv', sep='\t', index=False)
        zyg.to_csv(f'{args.out_prefix}_zygosity.tsv', sep='\t', index=False)
        descr.to_csv(f'{args.out_prefix}_descriptors.tsv', sep='\t', index=False)
        if thr is not None:
            thr.to_csv(f'{args.out_prefix}_thresholds.tsv', sep='\t', index=False)
        print(f'\nWrote {args.out_prefix}_{{calls,zygosity,descriptors'
              f'{",thresholds" if thr is not None else ""}}}.tsv')


if __name__ == '__main__':
    main()
