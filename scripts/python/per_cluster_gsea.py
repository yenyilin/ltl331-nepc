#!/usr/bin/env python3
"""
per_cluster_gsea.py — per-cluster Hallmark GSEA (NES) for Fig 3D, with selectable
gene RANKING because prerank GSEA depends entirely on the ranking, and different
rankings give different (both valid) PRAD-vs-NEPC patterns:

  --rank-mode vs_rest      cluster vs ALL other cells (Wilcoxon z). Measures
                           DISTINCTIVENESS: a program shared across all NEPC
                           clusters cancels within NEPC -> NEPC may not look warm.
  --rank-mode vs_reference cluster vs a fixed reference (default AR-high PRAD
                           11/15 = pre-castration origin; signed logFC of means).
                           Measures change-from-origin: NEPC lineage programs all
                           warm. Most likely to match a "per-cluster vs baseline"
                           Fig 3D.
  --rank-mode magnitude    rank by absolute mean expression in the cluster.
                           Measures activity magnitude: highly-expressed lineage
                           genes warm in all NEPC (housekeeping-biased).
  --rank-mode all          run all three (into <outdir>/<mode>/) so you can pick.

Why NES not module scores: rank-based, immune to the highly-expressed-gene bias.

Each mode writes per_cluster_NES.tsv (Fig 3D heatmap source) + per_cluster_FDR.tsv
+ kras_directional_NES.tsv. Feed the chosen NES to fig3d_nes_heatmap.py.

Deps: scanpy, anndata, gseapy, pandas, numpy, scipy.

Example:
  python per_cluster_gsea.py --h5ad ltl331_base.h5ad --cluster-key seurat_clusters \\
      --gmt h.all.v7.5.1.symbols.gmt --use-raw --rank-mode all \\
      --reference-clusters 11 15 --outdir per_cluster_gsea
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def cluster_rankings(adata, cluster_key, mode, use_raw, ref_mask):
    """Return {cluster_id: pd.Series(gene -> rank statistic)} for the given mode."""
    import scanpy as sc
    src = adata.raw if (use_raw and adata.raw is not None) else adata
    var = np.asarray(src.var_names)
    X = src.X
    clusters = adata.obs[cluster_key].astype(str).values
    uniq = sorted(set(clusters), key=lambda c: int(c) if c.isdigit() else 1e9)
    out = {}

    if mode == 'vs_rest':
        n_all = src.shape[1]
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon',
                                use_raw=use_raw, n_genes=n_all)
        rg = adata.uns['rank_genes_groups']
        for g in rg['names'].dtype.names:
            out[g] = pd.Series(np.asarray(rg['scores'][g], float),
                               index=np.asarray(rg['names'][g], str))
        return out

    def col_mean(mask):
        sub = X[mask]
        return np.asarray(sub.mean(axis=0)).ravel() if sp.issparse(sub) \
            else np.asarray(sub).mean(axis=0).ravel()

    if mode == 'vs_reference':
        mu_ref = col_mean(ref_mask)
        for c in uniq:
            m = clusters == c
            if m.sum() == 0:
                continue
            out[c] = pd.Series(col_mean(m) - mu_ref, index=var)   # signed logFC of lognorm
        return out

    if mode == 'magnitude':
        for c in uniq:
            m = clusters == c
            if m.sum() == 0:
                continue
            out[c] = pd.Series(col_mean(m), index=var)            # absolute activity
        return out

    raise ValueError(f'unknown rank-mode {mode}')


def run_gsea(rankings, gmt, up, dn, outdir, *, min_size, max_size, perms, threads, seed):
    import gseapy as gp
    outdir.mkdir(parents=True, exist_ok=True)
    nes_cols, fdr_cols, kras = {}, {}, []
    for g, rnk in rankings.items():
        rnk = rnk.dropna()
        rnk = rnk[~rnk.index.duplicated()].sort_values(ascending=False)
        pre = gp.prerank(rnk=rnk, gene_sets=gmt, min_size=min_size, max_size=max_size,
                         permutation_num=perms, threads=threads, seed=seed,
                         no_plot=True, outdir=None, verbose=False)
        res = pre.res2d.copy()
        res['Term'] = res['Term'].astype(str).str.replace(r'^.*__', '', regex=True)
        res = res.set_index('Term')
        nes = pd.to_numeric(res['NES'], errors='coerce')
        fcol = 'FDR q-val' if 'FDR q-val' in res.columns else 'fdr'
        fdr = pd.to_numeric(res[fcol], errors='coerce')
        nes_cols[g], fdr_cols[g] = nes, fdr
        u, d = float(nes.get(up, np.nan)), float(nes.get(dn, np.nan))
        kras.append({'cluster': g, 'NES_UP': u, 'FDR_UP': float(fdr.get(up, np.nan)),
                     'NES_DN': d, 'FDR_DN': float(fdr.get(dn, np.nan)),
                     'NES_UP_minus_DN': u - d})
        print(f'  cl {g:>2}: NES_UP={u:+.2f} NES_DN={d:+.2f} dir={u - d:+.2f}')
    pd.DataFrame(nes_cols).to_csv(outdir / 'per_cluster_NES.tsv', sep='\t')
    pd.DataFrame(fdr_cols).to_csv(outdir / 'per_cluster_FDR.tsv', sep='\t')
    pd.DataFrame(kras).sort_values('NES_UP_minus_DN', ascending=False).to_csv(
        outdir / 'kras_directional_NES.tsv', sep='\t', index=False)
    print(f'  -> wrote per_cluster_NES/FDR + kras_directional_NES to {outdir}/')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--gmt', required=True)
    ap.add_argument('--use-raw', dest='use_raw', action='store_true', default=True)
    ap.add_argument('--no-use-raw', dest='use_raw', action='store_false')
    ap.add_argument('--rank-mode', default='all',
                    choices=['vs_rest', 'vs_reference', 'magnitude', 'all'])
    ap.add_argument('--reference-clusters', nargs='+', default=['11', '15'],
                    help='reference cluster ids for vs_reference (default 11 15 = AR-high PRAD)')
    ap.add_argument('--reference-key', default=None,
                    help='alternative reference: obs column (e.g. timepoint)')
    ap.add_argument('--reference-value', default=None,
                    help='value in --reference-key defining the reference (e.g. preCx)')
    ap.add_argument('--up', default='HALLMARK_KRAS_SIGNALING_UP')
    ap.add_argument('--dn', default='HALLMARK_KRAS_SIGNALING_DN')
    ap.add_argument('--min-size', type=int, default=10)
    ap.add_argument('--max-size', type=int, default=500)
    ap.add_argument('--permutations', type=int, default=1000)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--outdir', default='per_cluster_gsea')
    args = ap.parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str).astype('category')

    # reference mask for vs_reference
    if args.reference_key and args.reference_value is not None:
        ref_mask = (adata.obs[args.reference_key].astype(str) == str(args.reference_value)).values
        ref_desc = f'{args.reference_key}={args.reference_value}'
    else:
        ref_mask = adata.obs[args.cluster_key].astype(str).isin(
            [str(c) for c in args.reference_clusters]).values
        ref_desc = f'clusters {args.reference_clusters}'

    modes = ['vs_rest', 'vs_reference', 'magnitude'] if args.rank_mode == 'all' \
        else [args.rank_mode]
    base = Path(args.outdir)
    for mode in modes:
        if mode == 'vs_reference':
            print(f'\n=== rank-mode: {mode} (reference = {ref_desc}, n={int(ref_mask.sum())}) ===')
        else:
            print(f'\n=== rank-mode: {mode} ===')
        rankings = cluster_rankings(adata, args.cluster_key, mode, args.use_raw, ref_mask)
        outdir = base / mode if len(modes) > 1 else base
        run_gsea(rankings, args.gmt, args.up, args.dn, outdir,
                 min_size=args.min_size, max_size=args.max_size,
                 perms=args.permutations, threads=args.threads, seed=args.seed)

    if len(modes) > 1:
        print(f'\nDone. Compare {base}/{{vs_rest,vs_reference,magnitude}}/kras_directional_NES.tsv '
              f'and pick the one matching the original Fig 3D; feed its per_cluster_NES.tsv '
              f'to fig3d_nes_heatmap.py.')


if __name__ == '__main__':
    main()
