#!/usr/bin/env python3
"""
terminal_reference_project.py — project the clinical cohorts ONTO the LTL331
reference manifold instead of forcing a symmetric co-embedding.

Why projection, not co-embedding: symmetric Harmony block-shifts each cohort and never
interleaves them (iLISI=1.0), and the 52k-cell PDX dominates the centroids (dispersion
collapses Gao/Li). Reference projection fixes both by construction because the LTL331 geometry
(which spans the full PRAD->bridge->NE trajectory incl. both NE terminals) is FIXED, and
clinical cells map in. The question becomes concrete: where on the LTL331 manifold do
patient NE cells land, and do they inherit the NE / ASCL1+/- terminal label?

Uses scanpy.tl.ingest (linear projection onto reference PCA + label transfer). lognorm
from .raw on the shared HVG set; no raw counts needed.

Outputs:
  confusion_<labelcol>.tsv      query original label x transferred label (the headline)
  query_ne_landing.tsv          for query NE cells: transferred-label breakdown per dataset
  projected_umap_coords.tsv     ref+query UMAP coords (ref umap, query ingested)
  projection_umap.{pdf,png}     ref grey + query colored by transferred label

Example:
  python scripts/terminal_reference_project.py --h5ad harmony.h5ad \
      --dataset-col dataset --ref-name LTL331 --label-col coarse_celltype \
      --split-ne --outdir data/terminal_projection
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def ne_subtype(raw, lab, ne_label):
    """Split NE into ASCL1+/- by ASCL1 lognorm; returns a modified label array."""
    out = lab.copy()
    if 'ASCL1' not in raw.var_names:
        print('[warn] ASCL1 not in .raw; skipping NE split')
        return out
    a = raw[:, 'ASCL1'].X
    a = np.asarray(a.todense()).ravel() if hasattr(a, 'todense') else np.asarray(a).ravel()
    ne = lab == ne_label
    out[ne] = np.where(a[ne] > 0, 'NE_ASCL1pos', 'NE_ASCL1neg')
    return out


def main():
    import anndata as ad
    import scanpy as sc

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--dataset-col', default='dataset')
    p.add_argument('--ref-name', default='LTL331')
    p.add_argument('--label-col', default='coarse_celltype')
    p.add_argument('--split-ne', action='store_true')
    p.add_argument('--ne-label', default='NE')
    p.add_argument('--hvg-flag', default='vst.variable')
    p.add_argument('--n-pcs', type=int, default=30)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default='data/terminal_projection')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    full = ad.read_h5ad(args.h5ad)
    if full.raw is None:
        raise SystemExit('[abort] need .raw lognorm; none found')

    # build a clean lognorm working object on the shared HVG set --------------
    work = full.raw.to_adata()
    work.obs = full.obs.copy()
    if args.hvg_flag in full.var.columns:
        hvg = [g for g in full.var_names[full.var[args.hvg_flag].astype(bool)]
               if g in work.var_names]
    else:
        hvg = list(full.var_names[full.var_names.isin(work.var_names)])
    work = work[:, hvg].copy()
    print(f'working object {work.shape} on {len(hvg)} HVGs (lognorm)')

    lab = work.obs[args.label_col].astype(str).values
    if args.split_ne:
        lab = ne_subtype(full.raw, lab, args.ne_label)
    work.obs['_label'] = lab
    ds = work.obs[args.dataset_col].astype(str).values

    ref = work[ds == args.ref_name].copy()
    qry = work[ds != args.ref_name].copy()
    if ref.n_obs == 0 or qry.n_obs == 0:
        raise SystemExit(f'[abort] ref/query empty (ref-name={args.ref_name!r}); '
                         f'datasets present: {sorted(set(ds))}')
    print(f'reference {args.ref_name}: {ref.n_obs} | query: {qry.n_obs} '
          f'({sorted(set(ds[ds!=args.ref_name]))})')

    # build the reference manifold (scale -> PCA -> neighbors -> UMAP) -------
    sc.pp.scale(ref, max_value=10)
    sc.tl.pca(ref, n_comps=args.n_pcs, random_state=args.seed)
    sc.pp.neighbors(ref, n_pcs=args.n_pcs, random_state=args.seed)
    sc.tl.umap(ref, random_state=args.seed)
    ref.obs['_label'] = ref.obs['_label'].astype('category')

    # project query onto the reference + transfer the label -----------------
    sc.tl.ingest(qry, ref, obs='_label')
    qry.obs['_transferred'] = qry.obs['_label_'] if '_label_' in qry.obs else qry.obs['_label']

    # (A) confusion: query ORIGINAL label x TRANSFERRED label ---------------
    conf = pd.crosstab(qry.obs[args.label_col].astype(str), qry.obs['_transferred'].astype(str))
    conf.to_csv(outdir / f'confusion_{args.label_col}.tsv', sep='\t')
    conf_frac = conf.div(conf.sum(1), axis=0).round(3)
    print(f'\n=== (A) query original {args.label_col} -> transferred label (row-fraction) ===')
    print(conf_frac.to_string())

    # (B) where do query NE cells land, per dataset -------------------------
    ne_orig = qry.obs[args.label_col].astype(str) == args.ne_label
    land = (qry.obs.loc[ne_orig.values]
            .assign(_ds=qry.obs.loc[ne_orig.values, args.dataset_col].astype(str),
                    _t=qry.obs.loc[ne_orig.values, '_transferred'].astype(str))
            .groupby(['_ds', '_t']).size().rename('n').reset_index())
    land.to_csv(outdir / 'query_ne_landing.tsv', sep='\t', index=False)
    print('\n=== (B) query NE cells: transferred-label landing per dataset ===')
    print(land.to_string(index=False))

    # coords + plot ---------------------------------------------------------
    coords = pd.concat([
        pd.DataFrame(ref.obsm['X_umap'], index=ref.obs_names, columns=['umap1', 'umap2'])
            .assign(role='reference', label=ref.obs['_label'].astype(str).values,
                    dataset=args.ref_name),
        pd.DataFrame(qry.obsm['X_umap'], index=qry.obs_names, columns=['umap1', 'umap2'])
            .assign(role='query', label=qry.obs['_transferred'].astype(str).values,
                    dataset=qry.obs[args.dataset_col].astype(str).values),
    ])
    coords.to_csv(outdir / 'projected_umap_coords.tsv', sep='\t')

    try:
        from _utils import setup_style, savefig
        plt = setup_style()
        fig, ax = plt.subplots(figsize=(5, 4.2))
        ax.scatter(coords.loc[coords.role == 'reference', 'umap1'],
                   coords.loc[coords.role == 'reference', 'umap2'],
                   s=2, c='0.82', linewidths=0, label='LTL331 ref')
        q = coords[coords.role == 'query']
        for lv in sorted(q.label.unique()):
            m = q.label == lv
            ax.scatter(q.loc[m, 'umap1'], q.loc[m, 'umap2'], s=4, linewidths=0, label=lv)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
        ax.legend(frameon=False, fontsize=6, markerscale=2, loc='best')
        ax.set_title('clinical cohorts projected onto LTL331')
        fig.tight_layout()
        savefig(fig, outdir, 'projection_umap', ('pdf', 'png'))
        plt.close(fig)
    except Exception as e:
        print(f'[plot] skipped: {e}')

    print(f'\nwrote outputs to {outdir}/')
    print('INTERPRET: terminal merge recovered if query NE cells transfer to NE '
          '(diagonal in A) and land on the LTL331 NE region (B), per cohort.')


if __name__ == '__main__':
    main()
