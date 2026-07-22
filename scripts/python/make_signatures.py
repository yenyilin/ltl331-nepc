#!/usr/bin/env python3
"""
make_signatures.py — build the LTL331 cluster signature file that feeds
patient_level_stats.py --signatures (Step 4 pseudo-bulk scoring).

Reads a Seurat FindAllMarkers / Supplementary Table 3 style DEG table and emits,
per cluster, the top-N up-regulated genes using the manuscript's own criteria
(log2FC > 0.5, FDR < 0.05, top 20) — the same signatures used to score the
clinical (Gao/Li) cells.

Input table needs columns for: gene, cluster, log2FC, adjusted p-value
(Seurat defaults: gene[/index], cluster, avg_log2FC, p_val_adj).

Output (matches patient_level_stats.py `_load_signatures`):
  TSV  -> one line per cluster:  <cluster><TAB>gene1,gene2,...,geneN
  JSON -> {"<cluster>": ["gene1", ...], ...}   (optional)

Deps: pandas only.

Example:
  python make_signatures.py \\
      --markers supp_table3_findallmarkers.csv \\
      --out ltl331_cluster_markers.tsv --out-json ltl331_cluster_markers.json \\
      --lfc-min 0.5 --padj-max 0.05 --top-n 20

Then:
  python patient_level_stats.py ... --signatures ltl331_cluster_markers.tsv
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--markers', required=True,
                   help='FindAllMarkers / Supp Table 3 (.csv/.tsv)')
    p.add_argument('--out', default='ltl331_cluster_markers.tsv')
    p.add_argument('--out-json', default=None)
    p.add_argument('--gene-col', default='gene')
    p.add_argument('--cluster-col', default='cluster')
    p.add_argument('--lfc-col', default='avg_log2FC')
    p.add_argument('--padj-col', default='p_val_adj')
    p.add_argument('--lfc-min', type=float, default=0.5)
    p.add_argument('--padj-max', type=float, default=0.05)
    p.add_argument('--top-n', type=int, default=20)
    p.add_argument('--sep', default=None, help='override delimiter (else by ext)')
    args = p.parse_args()

    # Sniff the delimiter from the header line so a tab-separated file with a
    # .csv extension (the Seurat FindAllMarkers export ships this way) parses
    # correctly without needing --sep on the command line.
    if args.sep:
        sep = args.sep
    else:
        with open(args.markers, 'r') as _fh:
            sep = '\t' if '\t' in _fh.readline() else ','
    df = pd.read_csv(args.markers, sep=sep)

    # gene may be a column or the index (Seurat often writes rownames)
    if args.gene_col not in df.columns:
        df = df.reset_index().rename(columns={'index': args.gene_col})
        if args.gene_col not in df.columns:
            raise SystemExit(f'gene column "{args.gene_col}" not found; '
                             f'columns: {list(df.columns)}')
    for col in (args.cluster_col, args.lfc_col, args.padj_col):
        if col not in df.columns:
            raise SystemExit(f'column "{col}" not found; columns: {list(df.columns)}')

    n0 = len(df)
    df = df[(df[args.lfc_col] > args.lfc_min) & (df[args.padj_col] < args.padj_max)]
    print(f'filtered {n0} -> {len(df)} rows (lfc>{args.lfc_min}, '
          f'FDR<{args.padj_max})')

    sigs = {}
    for clust, sub in df.groupby(args.cluster_col):
        top = (sub.sort_values(args.lfc_col, ascending=False)
                  .drop_duplicates(args.gene_col)
                  .head(args.top_n))
        genes = top[args.gene_col].astype(str).tolist()
        sigs[str(clust)] = genes

    # write TSV
    outp = Path(args.out)
    with open(outp, 'w') as fh:
        fh.write(f'# cluster <tab> comma-separated top genes '
                 f'(lfc>{args.lfc_min}, FDR<{args.padj_max}, top {args.top_n})\n')
        for clust in sorted(sigs, key=lambda c: int(c) if c.isdigit() else 1e9):
            fh.write(f'{clust}\t{",".join(sigs[clust])}\n')
    print(f'wrote {outp} ({len(sigs)} clusters)')

    if args.out_json:
        Path(args.out_json).write_text(json.dumps(sigs, indent=2))
        print(f'wrote {args.out_json}')

    # preview + warnings
    small = {c: len(g) for c, g in sigs.items() if len(g) < 3}
    if small:
        print(f'[warn] clusters with <3 genes (skipped in scoring): {small}')
    print('preview:')
    for clust in sorted(sigs, key=lambda c: int(c) if c.isdigit() else 1e9)[:5]:
        print(f'  {clust}: {sigs[clust][:8]}{" ..." if len(sigs[clust])>8 else ""}')


if __name__ == '__main__':
    main()
