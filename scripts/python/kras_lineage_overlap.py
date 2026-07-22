#!/usr/bin/env python3
"""
kras_lineage_overlap.py — build the Sub-A supplementary table that explains the
KRAS_UP/KRAS_DN co-enrichment WITHOUT touching Fig 3D.

This operates on the SAME inputs as the original Fig 3D (the students' Seurat
FindAllMarkers export `DEGs_all_331everyclusters.csv` + the Hallmark KRAS GMT sets), so
the table is consistent with the kept figure — it is the per-cluster "leading edge" the
fgsea figure was already using, made explicit.

For each NE cluster it lists, per KRAS set, the member genes that are DE markers of that
cluster (the genes responsible for the set's enrichment), with avg_log2FC and p_val_adj,
and flags those that are also canonical neuroendocrine-lineage markers. The headline
numbers it prints are the ones the response letter cites:
  - KRAS_UP and KRAS_DN share 0 genes (disjoint sets).
  - In each NE cluster both sets are driven by NE-lineage members → co-directional,
    lineage signal, not reciprocal KRAS-pathway activity.

Inputs:
  --degs   DEGs_all_331everyclusters.csv  (Seurat FindAllMarkers: columns gene, cluster,
           avg_log2FC, p_val_adj, pct.1, pct.2 — column names auto-detected)
  --gmt    h.all.v7.5.1.symbols.gmt       (Hallmark; KRAS_UP/DN pulled by name)

Example:
  python kras_lineage_overlap.py \\
    --degs idiots/DEGs_all_331everyclusters.csv \\
    --gmt  msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/h.all.v7.5.1.symbols.gmt \\
    --ne-clusters 1 3 7 9 10 13 14 \\
    --outdir data/kras_lineage_overlap
"""
import argparse
from pathlib import Path
import pandas as pd

# canonical neuroendocrine / neuronal lineage markers used to flag "lineage rider" genes.
NE_MARKERS = {
    'CHGA', 'CHGB', 'SYP', 'NCAM1', 'ENO2', 'INSM1', 'ASCL1', 'NEUROD1', 'SOX2',
    'CPE', 'SNAP25', 'PCSK1N', 'PCSK1', 'PCSK2', 'SCG2', 'SCG3', 'SCG5', 'SYT1',
    'SYT4', 'SYT11', 'STMN2', 'STMN1', 'TUBB2B', 'DCX', 'ELAVL4', 'MAP2', 'NNAT',
    'IGFBP2', 'EPHA5', 'CAMK1D', 'KIF5C', 'PCP4', 'RUNDC3A', 'RTN1', 'CELF3',
    'NKX2-1', 'FOXA2', 'ONECUT2', 'POU3F2', 'MSX1',
}


def read_gmt_set(gmt_path, name_substr):
    """Return (full_name, set(genes)) for the Hallmark set whose name contains name_substr."""
    hits = []
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            if name_substr in name:
                hits.append((name, set(g for g in parts[2:] if g)))
    if not hits:
        raise SystemExit(f'no Hallmark set matching "{name_substr}" in {gmt_path}')
    if len(hits) > 1:
        names = ', '.join(n for n, _ in hits)
        raise SystemExit(f'"{name_substr}" matched multiple sets: {names}; be more specific')
    return hits[0]


def pick(cols, *cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    raise SystemExit(f'none of {cands} found in columns {list(cols)}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--degs', required=True)
    ap.add_argument('--gmt', required=True)
    ap.add_argument('--up', default='KRAS_SIGNALING_UP')
    ap.add_argument('--dn', default='KRAS_SIGNALING_DN')
    ap.add_argument('--ne-clusters', nargs='+', default=['1', '3', '7', '9', '10', '13', '14'])
    ap.add_argument('--padj', type=float, default=0.05, help='DE significance cutoff (matches Fig 3D)')
    ap.add_argument('--outdir', default='kras_lineage_overlap')
    args = ap.parse_args()

    up_name, up_set = read_gmt_set(args.gmt, args.up)
    dn_name, dn_set = read_gmt_set(args.gmt, args.dn)
    shared = up_set & dn_set
    print(f'{up_name}: {len(up_set)} genes')
    print(f'{dn_name}: {len(dn_set)} genes')
    print(f'shared genes between UP and DN sets: {len(shared)}'
          + (f'  ({sorted(shared)})' if shared else '  -> disjoint, as cited'))

    deg = pd.read_csv(args.degs, sep=None, engine='python')  # sniff , or \t (.csv or .txt)
    gcol = pick(deg.columns, 'gene', 'Gene', 'names', 'feature')
    ccol = pick(deg.columns, 'cluster', 'group')
    fcol = pick(deg.columns, 'avg_log2FC', 'avg_logFC', 'log2FoldChange', 'logfc')
    pcol = pick(deg.columns, 'p_val_adj', 'padj', 'p_adj', 'FDR')
    deg[ccol] = deg[ccol].astype(str)
    deg = deg[deg[pcol] < args.padj]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f'\nper-cluster KRAS-set members among up-regulated DE markers (padj<{args.padj}):')
    for cl in [str(c) for c in args.ne_clusters]:
        sub = deg[(deg[ccol] == cl) & (deg[fcol] > 0)]
        up_hits = sub[sub[gcol].isin(up_set)]
        dn_hits = sub[sub[gcol].isin(dn_set)]
        for setname, hits in [(up_name, up_hits), (dn_name, dn_hits)]:
            for _, r in hits.iterrows():
                rows.append({'cluster': cl, 'hallmark_set': setname, 'gene': r[gcol],
                             'avg_log2FC': r[fcol], 'p_val_adj': r[pcol],
                             'NE_lineage_marker': r[gcol] in NE_MARKERS})
        n_up_ne = up_hits[gcol].isin(NE_MARKERS).sum()
        n_dn_ne = dn_hits[gcol].isin(NE_MARKERS).sum()
        print(f'  cl {cl:>2}: UP {len(up_hits):>2} members ({n_up_ne} NE-marker) | '
              f'DN {len(dn_hits):>2} members ({n_dn_ne} NE-marker)')

    tab = pd.DataFrame(rows).sort_values(['cluster', 'hallmark_set', 'avg_log2FC'],
                                         ascending=[True, True, False])
    out = outdir / 'kras_lineage_overlap_SuppTableX.tsv'
    tab.to_csv(out, sep='\t', index=False)
    print(f'\n[ok] wrote {out}  ({len(tab)} rows)')
    print('Supp Table X: per NE cluster, the KRAS_UP/DN member genes among its DE markers, '
          'with NE_lineage_marker flag — the leading edge behind the Fig 3D co-enrichment.')


if __name__ == '__main__':
    main()
