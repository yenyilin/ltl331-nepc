#!/usr/bin/env python3
"""
cluster17_transition_flux.py — NE bifurcation question, done the right way.

Quantify the direction of the cluster-17 bifurcation into the two NEPC branches,
ASCL1+ (cluster 10) and ASCL1- (via cluster 13 to cluster 7).
Absorption probabilities can't answer the ASCL1+ side cleanly because the ASCL1+ state
(cluster 10) is a TRANSIENT early-NE state, not a CellRank terminal — defining it as a
membership-terminal would create tautological fates (see methodology notes). Instead we
read the DIRECTED TRANSITION MATRIX and ask, for cells leaving cluster 17, how much
one-step probability mass flows to the ASCL1+ entry (cluster 10) versus the ASCL1- NE
entry (cluster 13 and the other NE clusters):

    flux[s -> t] = mean over cells i in cluster s of  sum_{j in cluster t} T[i, j]

This is a proper probability distribution over destination clusters per source cell,
averaged over the source cluster. We report:
  - raw outflux from cl17 to every cluster (rows sum to ~1; includes self & backward)
  - FORWARD-normalized branch preference among NE-entry clusters only (cl10 vs cl13/...)
    => "of the cells moving from cl17 into an early-NE state, what fraction head to the
        ASCL1+ entry vs the ASCL1- entry"

Transition matrix source, in priority order:
  1. --tmat-key  : a cell x cell matrix in adata.obsp (e.g. a saved CellRank kernel T,
                   or scVelo 'T_fwd'/velocity graph)
  2. --recompute : rebuild the documented kernel
                   (PseudotimeKernel(time=velocity_pt)*0.8 + ConnectivityKernel*0.2),
                   matching the absorption-probability Methods, and use its
                   .transition_matrix
  3. --paga      : fall back to directed PAGA cluster->cluster transitions in
                   adata.uns[<paga_key>]['transitions_confidence'] (already cluster-level)

Example (recompute, matches Methods):
  python cluster17_transition_flux.py --h5ad ltl331_velocity_cellrank.h5ad \\
      --cluster-key seurat_clusters --source 17 \\
      --ascl1-pos 10 --ascl1-neg-entry 13 --ne-clusters 1 3 7 9 13 14 \\
      --recompute --pt-key velocity_pt --conn-weight 0.2 \\
      --outdir data/cluster17_flux
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def cluster_flux(T, clusters, source, targets):
    """mean over source-cluster cells of summed transition mass into each target cluster."""
    src_rows = np.where(clusters == str(source))[0]
    if src_rows.size == 0:
        raise SystemExit(f'source cluster {source} has no cells')
    Tsrc = T[src_rows]                      # (n_src, n_cells)
    Tsrc = Tsrc.toarray() if sp.issparse(Tsrc) else np.asarray(Tsrc)
    rowsum = Tsrc.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    Tsrc = Tsrc / rowsum                    # guard against un-normalized rows
    out = {}
    for t in targets:
        cols = np.where(clusters == str(t))[0]
        out[str(t)] = float(Tsrc[:, cols].sum(axis=1).mean()) if cols.size else 0.0
    return out, src_rows.size


def get_tmat(adata, args):
    if args.tmat_key:
        if args.tmat_key not in adata.obsp:
            raise SystemExit(f'--tmat-key {args.tmat_key} not in adata.obsp ({list(adata.obsp)})')
        print(f'[tmat] using adata.obsp["{args.tmat_key}"]')
        return adata.obsp[args.tmat_key], 'cell'
    if args.recompute:
        import cellrank as cr
        print('[tmat] recomputing PseudotimeKernel*{:.1f} + ConnectivityKernel*{:.1f}'
              .format(1 - args.conn_weight, args.conn_weight))
        pk = cr.kernels.PseudotimeKernel(adata, time_key=args.pt_key).compute_transition_matrix()
        ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        comb = (1 - args.conn_weight) * pk + args.conn_weight * ck
        # a combined kernel already exposes .transition_matrix; compute only if absent
        if getattr(comb, 'transition_matrix', None) is None:
            comb.compute_transition_matrix()
        return comb.transition_matrix, 'cell'
    if args.paga:
        key = args.paga_key
        tc = adata.uns[key]['transitions_confidence']
        print(f'[tmat] using directed PAGA adata.uns["{key}"]["transitions_confidence"] (cluster-level)')
        return tc, 'cluster'
    raise SystemExit('provide one of --tmat-key / --recompute / --paga')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--source', default='17')
    ap.add_argument('--ascl1-pos', default='10', help='ASCL1+ early-NE entry cluster')
    ap.add_argument('--ascl1-neg-entry', default='13', help='ASCL1- NE entry cluster')
    ap.add_argument('--ne-clusters', nargs='+', default=['1', '3', '7', '9', '13', '14'],
                    help='all ASCL1- NE clusters (for the aggregate ASCL1- branch)')
    ap.add_argument('--tmat-key', default=None)
    ap.add_argument('--recompute', action='store_true')
    ap.add_argument('--pt-key', default='velocity_pt')
    ap.add_argument('--conn-weight', type=float, default=0.2)
    ap.add_argument('--paga', action='store_true')
    ap.add_argument('--paga-key', default='paga')
    ap.add_argument('--outdir', default='cluster17_flux')
    args = ap.parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    clusters = adata.obs[args.cluster_key].values
    all_clusters = sorted(set(clusters), key=lambda c: int(c) if c.isdigit() else 1e9)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    T, kind = get_tmat(adata, args)

    if kind == 'cluster':
        # PAGA: rows = source cluster index in adata.obs[cluster_key].cat.categories order
        cats = list(adata.obs[args.cluster_key].astype('category').cat.categories)
        Tc = T.toarray() if sp.issparse(T) else np.asarray(T)
        si = cats.index(str(args.source))
        row = Tc[si].copy()
        s = row.sum() or 1.0
        flux = {cats[j]: float(row[j] / s) for j in range(len(cats))}
        n_src = int((clusters == str(args.source)).sum())
    else:
        flux, n_src = cluster_flux(T, clusters, args.source, all_clusters)

    fl = pd.Series(flux).sort_values(ascending=False)
    fl.to_csv(outdir / f'cl{args.source}_outflux.tsv', sep='\t', header=['flux'])
    print(f'\ncluster {args.source} (n={n_src}) one-step out-flux to each cluster:')
    print(fl.round(4).to_string())

    pos = flux.get(str(args.ascl1_pos), 0.0)
    neg_entry = flux.get(str(args.ascl1_neg_entry), 0.0)
    neg_all = sum(flux.get(str(c), 0.0) for c in args.ne_clusters)
    # forward branch preference: ASCL1+ entry vs ASCL1- entry, renormalized to the two
    denom = pos + neg_entry
    pref_pos = pos / denom if denom else float('nan')
    print('\n=== bifurcation (transition flux out of cluster {}) ==='.format(args.source))
    print(f'  flux -> ASCL1+ entry (cl{args.ascl1_pos})        : {pos:.4f}')
    print(f'  flux -> ASCL1- NE entry (cl{args.ascl1_neg_entry}): {neg_entry:.4f}')
    print(f'  flux -> all ASCL1- NE {args.ne_clusters}: {neg_all:.4f}')
    print(f'  branch preference among the two entries: '
          f'ASCL1+ {pref_pos:.2f} vs ASCL1- {1 - pref_pos:.2f}'
          f'  (ratio ASCL1-/ASCL1+ = {neg_entry / pos:.1f}x)' if pos else '')
    pd.DataFrame([{'source': args.source, 'flux_ascl1pos_cl' + str(args.ascl1_pos): pos,
                   'flux_ascl1neg_entry_cl' + str(args.ascl1_neg_entry): neg_entry,
                   'flux_ascl1neg_all': neg_all, 'pref_ascl1pos': pref_pos,
                   'n_source_cells': n_src}]).to_csv(
        outdir / f'cl{args.source}_branch_preference.tsv', sep='\t', index=False)
    print(f'\n[ok] wrote outflux + branch preference to {outdir}/')


if __name__ == '__main__':
    main()
