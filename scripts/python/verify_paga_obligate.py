#!/usr/bin/env python3
"""
verify_paga_obligate.py — READ-ONLY check that cluster 17 is the obligate PAGA gateway.

Replicates the Fig 4A PAGA (paga_trajectory.py defaults: seurat_clusters, n_neighbors=8,
n_pcs=20) and tests, quantitatively:
  (1) the max PAGA connectivity between ANY adenocarcinoma cluster and ANY neuroendocrine
      cluster (cluster 17 excluded), or the "bypass" edge. Should be ~0 / below threshold.
  (2) cluster 17's connectivity to every cluster (its two bridges: a PRAD side + cluster 13).
  (3) cluster 13's connectivity to every cluster — its ONLY adenocarcinoma-side connector
      should be cluster 17.
  (4) cut-vertex test: with cluster 17 removed, are the PRAD and NE compartments disconnected
      at the figure threshold?

Writes nothing except stdout + an optional TSV. Usage:
  python verify_paga_obligate.py --h5ad /path/to/finaldong-log.h5ad [--threshold 0.09] [--out paga_connectivity.tsv]
"""
import argparse, json
import numpy as np

# compartments read off Fig 4A (the two node groups); cluster 17 is the bridge.
PRAD = {0, 2, 4, 5, 6, 8, 11, 12, 15, 16}
NE   = {1, 3, 7, 9, 10, 13, 14}
BRIDGE = 17

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--n-neighbors', type=int, default=8)
    ap.add_argument('--n-pcs', type=int, default=20)
    ap.add_argument('--threshold', type=float, default=0.09, help='Fig 4A drawn-edge threshold')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'])
    ap.add_argument('--recompute-pca', action='store_true')
    ap.add_argument('--out', default='paga_connectivity.tsv')
    args = ap.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.h5ad)
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        adata = adata[keep].copy()
    adata.obs['clusters'] = adata.obs[args.cluster_key].astype(str).astype('category')
    if args.recompute_pca or 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.paga(adata, groups='clusters')

    cats = list(adata.obs['clusters'].cat.categories)            # e.g. ['0','1',...,'17']
    C = adata.uns['paga']['connectivities'].toarray()            # symmetric, len(cats)^2
    idx = {int(c): i for i, c in enumerate(cats) if c.isdigit()}

    # full matrix -> TSV
    with open(args.out, 'w') as fh:
        fh.write('\t' + '\t'.join(cats) + '\n')
        for i, ci in enumerate(cats):
            fh.write(ci + '\t' + '\t'.join(f'{C[i, j]:.4f}' for j in range(len(cats))) + '\n')
    print(f'[wrote] {args.out}  ({len(cats)} clusters)\n')

    def conn(a, b):
        return C[idx[a], idx[b]] if a in idx and b in idx else float('nan')

    # (1) bypass edge: strongest PRAD<->NE connection, cluster 17 excluded
    bypass = [(p, n, conn(p, n)) for p in PRAD for n in NE if p in idx and n in idx]
    bypass = [x for x in bypass if not np.isnan(x[2])]
    bypass.sort(key=lambda x: -x[2])
    print('=== (1) Strongest direct PRAD<->NE edges (cluster 17 excluded) — the bypass test ===')
    for p, n, v in bypass[:6]:
        flag = '  <-- ABOVE THRESHOLD (bypass!)' if v >= args.threshold else ''
        print(f'   PRAD {p:>2} <-> NE {n:>2} : {v:.4f}{flag}')
    print(f'   => max PRAD-NE connectivity = {bypass[0][2]:.4f} (threshold {args.threshold}) '
          f'{"OBLIGATE OK: no bypass" if bypass[0][2] < args.threshold else "WARNING: a bypass edge exists"}\n')

    # (2) cluster 17 bridges
    print('=== (2) Cluster 17 connectivity (its bridges) ===')
    row17 = sorted(((int(c), conn(17, int(c))) for c in cats if c.isdigit() and int(c) != 17),
                   key=lambda x: -x[1])
    for c, v in row17[:8]:
        side = 'PRAD' if c in PRAD else ('NE' if c in NE else '?')
        print(f'   17 <-> {c:>2} ({side}): {v:.4f}')
    print()

    # (3) cluster 13 — its only PRAD-side connector should be 17
    print('=== (3) Cluster 13 connectivity — adenocarcinoma-side entry ===')
    prad_side = sorted(((c, conn(13, c)) for c in (PRAD | {17})), key=lambda x: -x[1])
    for c, v in prad_side[:6]:
        tag = 'cluster 17 (gateway)' if c == 17 else f'PRAD {c}'
        print(f'   13 <-> {tag}: {v:.4f}')
    top_prad = max(((c, conn(13, c)) for c in PRAD), key=lambda x: x[1])
    print(f'   => cl13 strongest *adenocarcinoma* connector = cluster {top_prad[0]} ({top_prad[1]:.4f}); '
          f'cl13<->17 = {conn(13,17):.4f}\n')

    # (4) cut-vertex test at threshold (remove 17, can we still reach NE from PRAD?)
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    keep_i = [i for i, c in enumerate(cats) if c.isdigit() and int(c) != 17]
    sub = (C[np.ix_(keep_i, keep_i)] >= args.threshold).astype(int)
    n_comp, labels = connected_components(sp.csr_matrix(sub), directed=False)
    lab = {int(cats[keep_i[k]]): labels[k] for k in range(len(keep_i))}
    prad_comp = {lab[c] for c in PRAD if c in lab}
    ne_comp = {lab[c] for c in NE if c in lab}
    print('=== (4) Cut-vertex test: remove cluster 17, threshold the graph ===')
    print(f'   PRAD clusters fall in component(s): {sorted(prad_comp)}')
    print(f'   NE   clusters fall in component(s): {sorted(ne_comp)}')
    disjoint = prad_comp.isdisjoint(ne_comp)
    print(f'   => {"CONFIRMED: removing cluster 17 DISCONNECTS PRAD from NE (cut vertex / obligate gateway)" if disjoint else "NOT a clean cut vertex — PRAD and NE share a component without 17"}')

if __name__ == '__main__':
    main()
