#!/usr/bin/env python3
"""
infercnv_clonal_check.py — exclude the "separate clonal origin" alternative (H3b)
for the cluster 17 -> ASCL1+ / ASCL1- bifurcation.

If cluster 17, cluster 10 (ASCL1+), and cluster 7 (ASCL1-) share the same
inferCNV-derived CNV architecture, they cannot be genetically distinct clones —
H3b is refuted. This script:

  1) reads an inferCNV-derived per-cell CNV matrix (cells x positions)
  2) joins it with cluster labels (from --h5ad or --labels-tsv)
  3) computes a cluster-cluster Pearson correlation matrix (means of per-cell
     CNV vectors) -> heatmap + TSV
  4) computes per-cell pairwise CNV correlations for {17-17, 17-10, 17-7, ...}
     -> distribution summaries + violins
  5) KS test: within-target distribution vs target-vs-target distribution
     (if indistinguishable -> shared clonal architecture)

Deps: numpy, pandas, scipy, matplotlib (anndata only if --h5ad used).

Example:
  python infercnv_clonal_check.py \\
      --cnv-matrix infercnv.observations.txt --transpose \\
      --h5ad ltl331.h5ad --cluster-key clusters \\
      --target-clusters 17 10 7 --baseline-clusters 11 15 \\
      --outdir infercnv_clonal --plot-format pdf png

Notes on the CNV matrix:
  - inferCNV typically exports its smoothed values as genes-as-rows x cells-as-cols.
    Pass --transpose to flip to cells-as-rows for this script.
  - Any per-cell CNV matrix works (e.g., HMM-discretised, expression-derived).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _utils import setup_style, savefig


def load_cnv(path, sep, transpose):
    s = sep or ('\t' if str(path).endswith(('.tsv', '.txt')) else ',')
    cnv = pd.read_csv(path, sep=s, index_col=0)
    if transpose:
        cnv = cnv.T
    print(f'CNV matrix: {cnv.shape[0]} cells x {cnv.shape[1]} positions')
    return cnv


def load_labels(args):
    if args.labels_tsv:
        df = pd.read_csv(args.labels_tsv, sep='\t', index_col=0)
        if args.cluster_key not in df.columns:
            # accept a single-column file
            if df.shape[1] == 1:
                df = df.rename(columns={df.columns[0]: args.cluster_key})
            else:
                raise SystemExit(
                    f'cluster column "{args.cluster_key}" not in labels-tsv; '
                    f'have: {list(df.columns)}')
        return df[args.cluster_key].astype(str)
    if args.h5ad:
        import anndata as ad
        a = ad.read_h5ad(args.h5ad)
        s = a.obs[args.cluster_key].astype(str)
        s.index = a.obs_names
        return s
    raise SystemExit('provide --labels-tsv or --h5ad for cluster labels')


def pairwise_pearson(mat):
    """Row-wise pairwise Pearson correlation of an (n x p) matrix, vectorised."""
    m = mat - mat.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(m, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    mn = m / norms
    return mn @ mn.T


def save_results_adata(args, outdir, *, sub_idx, lab_sub, cnv_sub, mean_cnv,
                       cluster_corr, cell_corr, pair_df, ks_rows):
    """Build a small results-only h5ad.

      X      (n_cells x n_clusters) = Pearson correlation of each subsampled
             cell's CNV vector to each cluster's mean CNV vector
      obs    cluster, within_cluster_mean_corr, top_cluster_match,
             top_cluster_corr
      var    n_cells_in_cluster, mean_abs_cnv (per cluster)
      uns    cluster_cnv_correlation, cell_pair_distribution_summary,
             ks_tests, infercnv_clonal_params

    The full input CNV matrix is intentionally NOT stored (kept compact).
    """
    try:
        import anndata as ad
    except Exception:
        print('[save] anndata not importable; skipping save')
        return

    A = cnv_sub.values.astype(float)
    B = mean_cnv.values.astype(float)
    Ac = A - A.mean(axis=1, keepdims=True)
    Bc = B - B.mean(axis=1, keepdims=True)
    An = Ac / (np.linalg.norm(Ac, axis=1, keepdims=True) + 1e-12)
    Bn = Bc / (np.linalg.norm(Bc, axis=1, keepdims=True) + 1e-12)
    cell_x_cluster = An @ Bn.T            # (n_cells x n_clusters)

    cluster_names = list(mean_cnv.index)

    # within-cluster mean correlation per cell (uses precomputed cell_corr)
    coh = np.full(len(sub_idx), np.nan)
    lab_arr = lab_sub.values
    cell_ids = np.asarray(sub_idx)
    for i, cid in enumerate(cell_ids):
        own = lab_sub.loc[cid]
        peers = cell_ids[(lab_arr == own) & (cell_ids != cid)]
        if len(peers) > 0:
            row = cell_corr.loc[cid, peers].values
            row = row[~np.isnan(row)]
            if len(row):
                coh[i] = float(row.mean())

    top_idx = np.nanargmax(cell_x_cluster, axis=1)
    obs = pd.DataFrame({
        'cluster': lab_sub.loc[cell_ids].values,
        'within_cluster_mean_corr': coh,
        'top_cluster_match': np.array(cluster_names)[top_idx],
        'top_cluster_corr': np.nanmax(cell_x_cluster, axis=1),
    }, index=cell_ids)

    var = pd.DataFrame({
        'n_cells_in_cluster': [int((lab_sub.values == c).sum())
                               for c in cluster_names],
        'mean_abs_cnv': [float(np.mean(np.abs(mean_cnv.loc[c].values)))
                         for c in cluster_names],
    }, index=cluster_names)

    adata_out = ad.AnnData(X=cell_x_cluster.astype('float32'), obs=obs, var=var)

    adata_out.uns['cluster_cnv_correlation'] = cluster_corr.to_dict()
    adata_out.uns['cell_pair_distribution_summary'] = pair_df.to_dict(orient='list')
    adata_out.uns['ks_tests'] = ks_rows
    adata_out.uns['infercnv_clonal_params'] = {
        'cnv_matrix': str(args.cnv_matrix),
        'transposed': bool(args.transpose),
        'cluster_source': (('h5ad: ' + args.h5ad) if args.h5ad
                           else ('labels: ' + str(args.labels_tsv))),
        'cluster_key': args.cluster_key,
        'target_clusters': [str(c) for c in args.target_clusters],
        'baseline_clusters': [str(c) for c in args.baseline_clusters],
        'max_per_cluster': int(args.max_per_cluster),
        'seed': int(args.seed),
        'note': ('Results-only AnnData. X = per-cell Pearson correlation to '
                 'each cluster-mean CNV vector. Full input CNV matrix not '
                 'stored.'),
    }

    suffix = args.save_suffix or 'clonal'
    out_path = Path(outdir) / f'infercnv_clonal_results_{suffix}.h5ad'
    adata_out.write_h5ad(out_path, compression='gzip')
    print(f'[save] wrote {out_path} '
          f'({adata_out.n_obs} cells x {adata_out.n_vars} clusters)')


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--cnv-matrix', required=True,
                   help='inferCNV cells x positions matrix (TSV/CSV); '
                        'pass --transpose if it is positions x cells')
    p.add_argument('--transpose', action='store_true')
    p.add_argument('--sep', default=None)
    p.add_argument('--labels-tsv', default=None,
                   help='cell_id<TAB>cluster   (alternative to --h5ad)')
    p.add_argument('--h5ad', default=None)
    p.add_argument('--cluster-key', default='clusters')
    p.add_argument('--target-clusters', nargs='+', default=['17', '10', '7'])
    p.add_argument('--baseline-clusters', nargs='+', default=['11', '15'],
                   help='AR-high PRAD clusters; in a monoclonal LTL331 these '
                        'should also correlate highly with the targets — '
                        'orthogonal sanity check')
    p.add_argument('--max-per-cluster', type=int, default=300)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--outdir', default='infercnv_clonal')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    p.add_argument('--save-adata', action='store_true',
                   help='write a small results-only h5ad: X = per-cell '
                        'correlation to each cluster mean CNV; full CNV matrix '
                        'not included')
    p.add_argument('--save-suffix', default=None,
                   help='filename suffix (default: clonal)')
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    cnv = load_cnv(args.cnv_matrix, args.sep, args.transpose)
    labels = load_labels(args)

    common = cnv.index.intersection(labels.index)
    if len(common) == 0:
        raise SystemExit('no overlapping cell IDs between CNV matrix and labels')
    cnv = cnv.loc[common]
    labels = labels.loc[common]
    print(f'aligned: {len(common)} cells')

    target = [str(c) for c in args.target_clusters]
    baseline = [str(c) for c in args.baseline_clusters]
    all_clusters = target + [c for c in baseline if c not in target]

    sub_idx = []
    for c in all_clusters:
        cells = labels.index[labels.values == c].tolist()
        if not cells:
            print(f'[warn] no cells in cluster {c}')
            continue
        n = min(args.max_per_cluster, len(cells))
        sub_idx.extend(list(rng.choice(cells, size=n, replace=False)))
    if not sub_idx:
        raise SystemExit('no cells selected; check cluster IDs')
    cnv_sub = cnv.loc[sub_idx]
    lab_sub = labels.loc[sub_idx]
    print(f'subsampled to {len(sub_idx)} cells across '
          f'{lab_sub.nunique()} clusters '
          f'(max {args.max_per_cluster}/cluster)')

    # ---- (1) cluster-level: mean CNV per cluster, cluster-cluster Pearson ----
    mean_cnv = cnv_sub.groupby(lab_sub.values).mean()
    cluster_corr = mean_cnv.T.corr()
    cluster_corr.to_csv(outdir / 'cluster_cnv_correlation.tsv', sep='\t')
    print('\ncluster-mean CNV Pearson correlation:')
    print(cluster_corr.round(3).to_string())

    # ---- (2) per-cell pairwise Pearson on the subsampled matrix ----
    mat = cnv_sub.values.astype(float)
    cc = pairwise_pearson(mat)
    np.fill_diagonal(cc, np.nan)
    cell_corr = pd.DataFrame(cc, index=sub_idx, columns=sub_idx)

    cells_by = {c: lab_sub.index[lab_sub.values == c].tolist()
                for c in all_clusters if (lab_sub.values == c).any()}
    pair_dists = {}
    rows = []
    for c1 in target:
        if c1 not in cells_by:
            continue
        for c2 in all_clusters:
            if c2 not in cells_by:
                continue
            block = cell_corr.loc[cells_by[c1], cells_by[c2]].values
            block = block[~np.isnan(block)]
            pair_dists[(c1, c2)] = block
            if len(block):
                rows.append({'A': c1, 'B': c2, 'n_pairs': len(block),
                             'mean_r': float(np.mean(block)),
                             'median_r': float(np.median(block)),
                             'sd_r': float(np.std(block))})
    pair_df = pd.DataFrame(rows)
    pair_df.to_csv(outdir / 'cell_cell_correlation_distributions.tsv',
                   sep='\t', index=False)
    print('\nper-cluster-pair cell-cell correlation summary:')
    print(pair_df.round(3).to_string(index=False))

    # ---- (3) KS test: within-target baseline vs target-target ----
    print('\nKS test — within-target baseline vs target-vs-other-target:')
    ks_rows = []
    for c1 in target:
        within = pair_dists.get((c1, c1))
        if within is None or len(within) < 5:
            continue
        for c2 in target:
            if c1 == c2:
                continue
            between = pair_dists.get((c1, c2))
            if between is None or len(between) < 5:
                continue
            ks, pv = stats.ks_2samp(within, between)
            verdict = ('indistinguishable (clones merge)'
                       if ks < 0.1 else 'distinguishable (note effect size)')
            print(f'   {c1}-{c1} vs {c1}-{c2}: KS={ks:.3f} p={pv:.3g}  -> {verdict}')
            ks_rows.append({'within': c1, 'between': c2,
                            'KS': float(ks), 'p': float(pv),
                            'verdict': verdict})

    # ---- (4) plots ----
    plt = setup_style()
    n = len(cluster_corr)
    fig, ax = plt.subplots(figsize=(0.55 * n + 1.5, 0.55 * n + 1.5))
    im = ax.imshow(cluster_corr.values, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(cluster_corr.columns, rotation=90, fontsize=10)
    ax.set_yticks(range(n)); ax.set_yticklabels(cluster_corr.index, fontsize=10)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f'{cluster_corr.values[i, j]:.2f}',
                    ha='center', va='center', fontsize=9,
                    color='white' if abs(cluster_corr.values[i, j]) > 0.5 else 'black')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label('Pearson r', fontsize=9)
    ax.set_title('cluster mean CNV correlation', fontsize=11)
    savefig(fig, outdir, 'cluster_cnv_correlation_heatmap', args.plot_format)
    plt.close(fig)

    inter = target[0]
    if inter in cells_by:
        labs, data = [], []
        for c2 in all_clusters:
            d = pair_dists.get((inter, c2))
            if d is not None and len(d):
                labs.append(f'{inter}-{c2}')
                data.append(d)
        if labs:
            fig, ax = plt.subplots(figsize=(max(5, 1.1 * len(labs) + 2), 3.5))
            ax.violinplot(data, showmeans=True)
            ax.set_xticks(range(1, len(labs) + 1))
            ax.set_xticklabels(labs, fontsize=7)
            ax.set_ylabel('cell-cell CNV Pearson r')
            ax.set_title(f'Per-cell CNV correlation distributions vs cluster {inter}')
            savefig(fig, outdir, f'cnv_correlation_vs_{inter}', args.plot_format)
            plt.close(fig)

    # ---- (4b) optional save of a small results-only h5ad ----
    if args.save_adata:
        print('\n[save] writing results-only h5ad ...')
        save_results_adata(args, outdir,
                           sub_idx=sub_idx, lab_sub=lab_sub, cnv_sub=cnv_sub,
                           mean_cnv=mean_cnv, cluster_corr=cluster_corr,
                           cell_corr=cell_corr, pair_df=pair_df,
                           ks_rows=ks_rows)

    # ---- (5) interpretation ----
    print('\nINTERPRET (H3b "separate clonal origin"):')
    print('  - cluster-mean correlations near 1 across targets => same clone.')
    print('  - within-target (e.g. 17-17) and target-vs-target (e.g. 17-7, 17-10)')
    print('    per-cell distributions overlapping (KS small) => clones not')
    print('    distinguishable by inferCNV => H3b NOT supported.')
    print('  - PRAD baseline clusters correlating highly with the NEPC targets')
    print('    is consistent with the manuscript inferCNV claim (conserved CNV')
    print('    across all clusters), independently refuting H3b.')
    print('  - Definitive clonal exclusion would require WES/scWGS at single-cell')
    print('    resolution; inferCNV is a strong but expression-derived proxy.')
    print('\nDone.')


if __name__ == '__main__':
    main()
