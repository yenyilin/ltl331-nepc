#!/usr/bin/env python3
"""
terminal_metaneighbor.py — cross-cohort replicability of the TERMINAL states
without trusting the failed co-embedding (iLISI=1.0).

MetaNeighbor-style unsupervised neighbor-voting AUROC: for every (dataset | cell_type)
group, train neighbor-voting on that group and test how well it recovers the SAME cell
type in the OTHER datasets (one-vs-rest within each test dataset). High cross-dataset
same-type AUROC = the state replicates across cohorts. This is embedding-free, so it is
immune to the iLISI=1.0 problem. It answers "is LTL331 NE recognizable as Gao/Li NE?"

Primary label = coarse_celltype (PRAD/NE/Basal/Other, already in obs). Optionally also
splits NE into ASCL1+/ASCL1- (by ASCL1 lognorm) to test whether the terminal NE SUBtypes
replicate, not just "NE".

Reads lognorm from .raw (29,622 genes); scores on the shared HVG set. Subsamples per
group to keep the cell-cell network tractable.

Example:
  python scripts/terminal_metaneighbor.py --h5ad harmony.h5ad \
      --dataset-col dataset --label-col coarse_celltype \
      --split-ne --max-per-group 300 --outdir data/terminal_metaneighbor
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def auroc(scores, pos_mask):
    from scipy.stats import rankdata
    pos_mask = np.asarray(pos_mask, bool)
    n_pos, n_neg = int(pos_mask.sum()), int((~pos_mask).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    r = rankdata(scores)
    return (r[pos_mask].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def rank_columns(M):
    """Per-column rank in [0,1] (MetaNeighbor ranks the cell-cell network)."""
    order = np.argsort(M, axis=0)
    ranks = np.empty_like(M)
    ar = np.arange(M.shape[0])
    for k in range(M.shape[1]):
        ranks[order[:, k], k] = ar
    return ranks / max(M.shape[0] - 1, 1)


def main():
    import anndata as ad

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--dataset-col', default='dataset')
    p.add_argument('--label-col', default='coarse_celltype')
    p.add_argument('--split-ne', action='store_true',
                   help='also split NE into ASCL1+/ASCL1- by ASCL1 lognorm')
    p.add_argument('--ne-label', default='NE')
    p.add_argument('--hvg-flag', default='vst.variable',
                   help='var column flagging HVGs; falls back to all shared genes')
    p.add_argument('--max-per-group', type=int, default=300)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default='data/terminal_metaneighbor')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    adata = ad.read_h5ad(args.h5ad)
    if adata.raw is None:
        raise SystemExit('[abort] need .raw lognorm; none found')
    raw = adata.raw
    obs = adata.obs

    # labels --------------------------------------------------------------
    ds = obs[args.dataset_col].astype(str).values
    lab = obs[args.label_col].astype(str).values.copy()
    if args.split_ne:
        if 'ASCL1' not in raw.var_names:
            print('[warn] ASCL1 not in .raw var; skipping NE split')
        else:
            ascl1 = raw[:, 'ASCL1'].X
            ascl1 = np.asarray(ascl1.todense()).ravel() if hasattr(ascl1, 'todense') \
                else np.asarray(ascl1).ravel()
            ne = lab == args.ne_label
            lab[ne] = np.where(ascl1[ne] > 0, 'NE_ASCL1pos', 'NE_ASCL1neg')
            print(f'[ne-split] NE -> ASCL1pos={int((lab=="NE_ASCL1pos").sum())} '
                  f'ASCL1neg={int((lab=="NE_ASCL1neg").sum())}')

    # subsample per (dataset|label) group --------------------------------
    groups = pd.Series([f'{d}|{l}' for d, l in zip(ds, lab)])
    keep_idx = []
    for g, idx in groups.groupby(groups).groups.items():
        idx = np.asarray(idx)
        if len(idx) > args.max_per_group:
            idx = rng.choice(idx, args.max_per_group, replace=False)
        keep_idx.append(idx)
    keep = np.sort(np.concatenate(keep_idx))
    print(f'subsampled {len(keep)} cells across {groups.nunique()} (dataset|label) groups')

    # variable-gene lognorm matrix ---------------------------------------
    if args.hvg_flag in adata.var.columns:
        hvg = adata.var_names[adata.var[args.hvg_flag].astype(bool)]
        genes = [g for g in hvg if g in raw.var_names]
    else:
        genes = list(adata.var_names[adata.var_names.isin(raw.var_names)])
    print(f'using {len(genes)} shared variable genes')
    gi = [raw.var_names.get_loc(g) for g in genes]
    X = raw.X[keep][:, gi]
    X = np.asarray(X.todense()) if hasattr(X, 'todense') else np.asarray(X)
    sub_ds, sub_lab = ds[keep], lab[keep]
    sub_grp = np.array([f'{d}|{l}' for d, l in zip(sub_ds, sub_lab)])

    # Spearman cell-cell network: rank genes across cells, zscore, dot -----
    from scipy.stats import rankdata
    Xr = np.apply_along_axis(rankdata, 0, X)
    Xr = (Xr - Xr.mean(0)) / (Xr.std(0) + 1e-9)
    net = (Xr @ Xr.T) / Xr.shape[1]            # [-1,1] Spearman
    net = rank_columns(net)                     # MetaNeighbor column-rank -> [0,1]
    deg = net.sum(1, keepdims=True) + 1e-9

    # neighbor voting: AUROC[train group i, test group j] (j's study, one-vs-rest)
    grp_names = sorted(set(sub_grp))
    Lmat = np.stack([(sub_grp == g).astype(float) for g in grp_names], axis=1)
    votes = (net @ Lmat) / deg                  # cells x train-groups
    auroc_mat = pd.DataFrame(np.nan, index=grp_names, columns=grp_names)
    for j, gj in enumerate(grp_names):
        ds_j = gj.split('|')[0]
        test = sub_ds == ds_j                    # one-vs-rest WITHIN test study
        pos = (sub_grp == gj) & test
        for i, gi_ in enumerate(grp_names):
            if gi_.split('|')[0] == ds_j:        # train must be a DIFFERENT study
                continue
            auroc_mat.loc[gi_, gj] = auroc(votes[test, i], pos[test])
    auroc_mat.to_csv(outdir / 'metaneighbor_auroc_matrix.tsv', sep='\t')

    # cross-dataset same-type summary ------------------------------------
    rows = []
    cell_types = sorted({g.split('|')[1] for g in grp_names})
    for ct in cell_types:
        vals = []
        for gi_ in grp_names:
            for gj in grp_names:
                if (gi_.split('|')[1] == ct and gj.split('|')[1] == ct
                        and gi_.split('|')[0] != gj.split('|')[0]):
                    v = auroc_mat.loc[gi_, gj]
                    if not np.isnan(v):
                        vals.append(v)
        rows.append({'cell_type': ct, 'mean_cross_dataset_AUROC': np.mean(vals) if vals else np.nan,
                     'min': np.min(vals) if vals else np.nan,
                     'max': np.max(vals) if vals else np.nan, 'n_pairs': len(vals)})
    summ = pd.DataFrame(rows).sort_values('mean_cross_dataset_AUROC', ascending=False)
    summ.to_csv(outdir / 'metaneighbor_crossdataset_summary.tsv', sep='\t', index=False)

    print('\n=== cross-dataset same-type replicability (neighbor-voting AUROC) ===')
    print(summ.round(3).to_string(index=False))
    print('\n  AUROC > 0.80 = the state is replicable across cohorts; ~0.5 = not.')
    print(f'\nwrote matrix + summary to {outdir}/')


if __name__ == '__main__':
    main()
