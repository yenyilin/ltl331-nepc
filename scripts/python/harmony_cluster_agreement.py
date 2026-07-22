#!/usr/bin/env python3
"""
harmony_cluster_agreement.py

Compare post-Harmony clustering against per-dataset `original_cluster` labels
for LTL331 + Gao + Li integration.

For each Leiden resolution requested:
  - Computes ARI / NMI / homogeneity / completeness per cohort (Gao, Li, LTL331)
    plus an ALL row across all cells
  - Writes per-cohort confusion matrices
  - Extracts the cluster-17 absorption breakdown (which new clusters absorbed
    the rare transient cluster 17 cells from LTL331)

Outputs in OUTDIR:
  metrics_summary.tsv             long-format master (resolution × group × metrics)
  metrics_wide.tsv                wide pivot for direct response-letter table
  metrics_r<R>.tsv                per-resolution snapshot
  confusion_<group>_r<R>.tsv      per-cohort confusion matrix
  cluster<X>_absorption.tsv       cluster-17 (or other target) absorption across resolutions

Usage:
    python harmony_cluster_agreement.py \\
        --h5ad harmony.h5ad \\
        --outdir ./agreement \\
        --resolutions 0.4 0.5 0.55 0.6 0.65 0.8 \\
        --group-col group \\
        --original-col original_cluster \\
        --target-cluster 17
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
)

from _utils import setup_style, savefig as _save_fig


def plot_metrics_sweep(master, outdir, formats):
    """Line plot of ARI/NMI/homogeneity/completeness vs resolution per cohort."""
    plt = setup_style()

    metrics = ['ARI', 'NMI', 'homogeneity', 'completeness']
    groups = [g for g in master['group'].unique() if g != 'ALL'] + ['ALL']

    fig, axes = plt.subplots(1, 4, figsize=(11, 2.6), sharex=True)
    for ax, metric in zip(axes, metrics):
        for grp in groups:
            sub = master[master['group'] == grp].sort_values('resolution')
            ax.plot(sub['resolution'], sub[metric],
                    marker='o', markersize=3, linewidth=1,
                    label=grp,
                    linestyle='--' if grp == 'ALL' else '-')
        ax.set_title(metric)
        ax.set_xlabel('Leiden resolution')
        ax.grid(True, alpha=0.3, linewidth=0.4)
    axes[0].set_ylabel('score')
    axes[-1].legend(frameon=False, fontsize=6, loc='best')
    fig.tight_layout()

    _save_fig(fig, outdir, 'metrics_sweep', formats)
    plt.close(fig)


def plot_target_absorption(c_df, outdir, formats, target_cluster):
    """Heatmap of how target-cluster cells distribute across new clusters.

    c_df: rows = new Leiden clusters, cols = resolutions, values = cell counts.
    Annotates counts; columns are column-normalised for the colour scale so
    the dominant absorbing clusters are visible regardless of total n.
    """
    plt = setup_style()

    # Drop new clusters that never absorb any target cells (all-zero rows)
    c_df = c_df.loc[(c_df.sum(axis=1) > 0)]
    if c_df.empty:
        print('  [plot] no target-cluster absorption to plot')
        return

    col_frac = c_df.div(c_df.sum(axis=0).replace(0, np.nan), axis=1)

    fig, ax = plt.subplots(
        figsize=(0.7 * c_df.shape[1] + 1.5, 0.32 * c_df.shape[0] + 1.5)
    )
    im = ax.imshow(col_frac.values, aspect='auto', cmap='Reds',
                   vmin=0, vmax=1)

    ax.set_xticks(range(c_df.shape[1]))
    ax.set_xticklabels(c_df.columns, rotation=0)
    ax.set_yticks(range(c_df.shape[0]))
    ax.set_yticklabels(c_df.index)
    ax.set_xlabel('Leiden resolution')
    ax.set_ylabel('new cluster')
    ax.set_title(f'Original cluster {target_cluster} absorption')

    # Annotate raw counts
    for i in range(c_df.shape[0]):
        for j in range(c_df.shape[1]):
            n = c_df.values[i, j]
            if n > 0:
                frac = col_frac.values[i, j]
                ax.text(j, i, str(int(n)), ha='center', va='center',
                        fontsize=6,
                        color='white' if frac > 0.5 else 'black')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('fraction of target cells (per resolution)', fontsize=6)
    fig.tight_layout()

    _save_fig(fig, outdir, f'cluster{target_cluster}_absorption', formats)
    plt.close(fig)


def compute_or_load_leiden(adata, resolutions, key_prefix='leiden_r',
                            cluster_key='harmony', n_neighbors=30,
                            random_state=0):
    """Ensure a Leiden column exists for each resolution; compute if missing.

    Reuses an existing neighbor graph ONLY if its connectivities matrix is actually
    present in `adata.obsp` (a populated `uns['neighbors']` whose `.obsp` graph was
    dropped during h5ad round-trip is NOT usable — Leiden would crash in _choose_graph).
    Otherwise (re)builds on the Harmony embedding (`X_harmony_dataset` preferred), and
    LOUDLY warns if it can only fall back to PCA (= naive, un-integrated space).
    """
    nb = adata.uns.get('neighbors', {})
    ck = nb.get('connectivities_key', 'connectivities') if isinstance(nb, dict) else 'connectivities'
    if 'neighbors' in adata.uns and ck in adata.obsp:
        print(f'[neighbors] using existing graph (obsp[{ck!r}], '
              f'{adata.obsp[ck].shape[0]:,} cells)')
    else:
        if 'neighbors' in adata.uns:
            print(f'[neighbors] uns["neighbors"] present but obsp[{ck!r}] missing '
                  f'(dropped on save) — rebuilding')
        # prefer the integrated Harmony rep; PCA is naive and only a last resort
        candidates = ['X_harmony_dataset', f'X_{cluster_key}', cluster_key, 'X_harmony', 'X_pca']
        rep = next((c for c in candidates if c in adata.obsm), None)
        if rep is None:
            raise ValueError(f'No embedding found in adata.obsm. Tried: {candidates}; '
                             f'have {list(adata.obsm)}')
        if rep == 'X_pca':
            print('[neighbors] *** WARN: no Harmony embedding found — building on X_pca '
                  '(NAIVE/un-integrated space; clusters will NOT reflect integration) ***')
        print(f'[neighbors] building graph on {rep} (k={n_neighbors})')
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=n_neighbors)

    for r in resolutions:
        key = f'{key_prefix}{r}'
        if key not in adata.obs.columns:
            print(f'[leiden] computing resolution {r} -> {key} '
                  f'(random_state={random_state})')
            sc.tl.leiden(adata, resolution=r, key_added=key,
                         random_state=random_state)
        else:
            print(f'[leiden] reusing existing column {key}')


def cluster_agreement(adata, group_col, a, b):
    """Per-cohort + ALL agreement metrics + per-cohort confusion matrices.

    Coerces all label columns to str to avoid categorical-dtype assignment
    errors and ensures sklearn metrics get hashable labels.
    """
    df = adata.obs[[group_col, a, b]].dropna().copy()
    df[group_col] = df[group_col].astype(str)
    df[a] = df[a].astype(str)
    df[b] = df[b].astype(str)

    rows = []
    confusions = {}

    for grp, sub in df.groupby(group_col):
        n = len(sub)
        if n < 2:
            rows.append({'group': grp, 'n': n,
                         'ARI': np.nan, 'NMI': np.nan,
                         'homogeneity': np.nan, 'completeness': np.nan})
            continue
        la, lb = sub[a], sub[b]
        rows.append({
            'group': grp, 'n': n,
            'ARI':          adjusted_rand_score(la, lb),
            'NMI':          normalized_mutual_info_score(la, lb),
            'homogeneity':  homogeneity_score(la, lb),
            'completeness': completeness_score(la, lb),
        })
        confusions[grp] = pd.crosstab(la, lb, margins=True,
                                       rownames=[a], colnames=[b])

    la_all, lb_all = df[a], df[b]
    rows.append({
        'group': 'ALL', 'n': len(df),
        'ARI':          adjusted_rand_score(la_all, lb_all),
        'NMI':          normalized_mutual_info_score(la_all, lb_all),
        'homogeneity':  homogeneity_score(la_all, lb_all),
        'completeness': completeness_score(la_all, lb_all),
    })

    return pd.DataFrame(rows).set_index('group'), confusions


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument('--h5ad', required=True,
                   help='Input AnnData (post-Harmony)')
    p.add_argument('--outdir', default='./agreement',
                   help='Output directory for TSVs')
    p.add_argument('--resolutions', nargs='+', type=float,
                   default=[0.4, 0.5, 0.6, 0.8],
                   help='Leiden resolutions to sweep')
    p.add_argument('--group-col', default='group',
                   help='obs column with dataset labels (e.g., Gao/Li/LTL331)')
    p.add_argument('--original-col', default='original_cluster',
                   help='obs column with per-dataset pre-Harmony cluster IDs')
    p.add_argument('--cluster-key', default='harmony',
                   help='Embedding key in adata.obsm to use for neighbors')
    p.add_argument('--leiden-prefix', default='leiden_r',
                   help='Prefix for new Leiden columns')
    p.add_argument('--random-state', type=int, default=0,
                   help='Seed for sc.tl.leiden (Leiden is stochastic; pin for '
                        'reproducible ARI/NMI sweeps)')
    p.add_argument('--target-cluster', default='17',
                   help='Cluster ID in --original-col to extract '
                        'absorption breakdown for (default: 17)')
    p.add_argument('--ltl331-label', default='LTL331',
                   help='Value in --group-col that identifies LTL331 cells')
    p.add_argument('--plots', action='store_true',
                   help='Also render figures (metrics sweep + target absorption)')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'],
                   help='Figure format(s) to write when --plots is set')
    p.add_argument('--save-clusters', default=None, metavar='PATH',
                   help='Write a SLIM obs-only h5ad (all obs incl. the clustering '
                        'columns; no X/raw/obsm — a few MB) so compare_resolution_'
                        'trajectory.py can consume it as --h5ad. Avoids rewriting the '
                        'full multi-GB object.')
    p.add_argument('--use-existing-cols', action='store_true',
                   help='Score PRE-EXISTING clustering columns instead of computing Leiden '
                        '(e.g. the Seurat/Louvain RNA_snn_res.* columns, for consistency '
                        'with the manuscript clustering). Column for resolution r is '
                        '"<existing-prefix><r>".')
    p.add_argument('--existing-prefix', default='RNA_snn_res.',
                   help='Prefix for pre-existing cluster columns (with --use-existing-cols); '
                        'default "RNA_snn_res." (Seurat). Column = prefix + resolution.')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    print(f'Loaded: {adata.n_obs} cells x {adata.n_vars} genes')
    print(f'  group counts: {dict(adata.obs[args.group_col].value_counts())}')

    if args.original_col not in adata.obs.columns:
        raise ValueError(
            f'--original-col "{args.original_col}" not in adata.obs '
            f'(columns: {list(adata.obs.columns)})'
        )

    if args.use_existing_cols:
        prefix = args.existing_prefix
        missing = [r for r in args.resolutions if f'{prefix}{r}' not in adata.obs.columns]
        if missing:
            raise SystemExit(
                f'[abort] --use-existing-cols: columns {[f"{prefix}{r}" for r in missing]} '
                f'not in obs. Present {prefix}* cols: '
                f'{[c for c in adata.obs.columns if str(c).startswith(prefix)]}')
        print(f'[clusters] scoring PRE-EXISTING columns {prefix}{{res}} (no Leiden computed)')
    else:
        compute_or_load_leiden(
            adata, args.resolutions,
            key_prefix=args.leiden_prefix,
            cluster_key=args.cluster_key,
            random_state=args.random_state,
        )
        prefix = args.leiden_prefix

    if args.save_clusters:
        import anndata as ad
        cluster_cols = [f'{prefix}{r}' for r in args.resolutions]
        slim = ad.AnnData(obs=adata.obs.copy())   # obs-only: no X/raw/obsm -> small + portable
        slim.write_h5ad(args.save_clusters)
        print(f'[save] slim obs-only object -> {args.save_clusters} '
              f'({slim.n_obs:,} cells; cluster cols: {cluster_cols}) '
              f'-> feed to compare_resolution_trajectory.py --h5ad')

    long_rows = []
    target_breakdown = {}

    for r in args.resolutions:
        leiden_key = f'{prefix}{r}'
        n_clusters = adata.obs[leiden_key].nunique()

        print(f'\n=== resolution {r}  (n_clusters = {n_clusters}) ===')
        metrics, confusions = cluster_agreement(
            adata, args.group_col, args.original_col, leiden_key
        )
        metrics.insert(0, 'resolution', r)
        metrics.insert(1, 'n_clusters', n_clusters)
        print(metrics.round(4).to_string())

        # Per-cohort confusion matrices
        for grp, cm in confusions.items():
            cm.to_csv(outdir / f'confusion_{grp}_r{r}.tsv', sep='\t')

        # Per-resolution metrics snapshot
        metrics.reset_index().to_csv(
            outdir / f'metrics_r{r}.tsv', sep='\t', index=False,
        )

        # Append for long-format master
        for grp, row in metrics.iterrows():
            long_rows.append({
                'resolution': r,
                'n_clusters': n_clusters,
                'group': grp,
                **{k: row[k] for k in ['n', 'ARI', 'NMI',
                                        'homogeneity', 'completeness']},
            })

        # Target-cluster absorption (LTL331 cells only)
        ltl_mask = adata.obs[args.group_col].astype(str) == args.ltl331_label
        target_mask = ltl_mask & (
            adata.obs[args.original_col].astype(str) == args.target_cluster
        )
        n_target = int(target_mask.sum())
        if n_target > 0:
            breakdown = adata.obs.loc[target_mask, leiden_key].value_counts()
            target_breakdown[r] = breakdown
            print(f'\nCluster {args.target_cluster} (n={n_target}) '
                  f'absorption at r={r}:')
            print(breakdown.to_string())
        else:
            print(f'\n[warn] no cells with {args.group_col}=={args.ltl331_label} '
                  f'and {args.original_col}=={args.target_cluster}')

    # Long-format master across all resolutions
    master = pd.DataFrame(long_rows)
    master.to_csv(outdir / 'metrics_summary.tsv', sep='\t', index=False)
    print(f'\nSaved long-format master: {outdir/"metrics_summary.tsv"}')

    # Wide-format pivot for direct response-letter table
    wide = master.pivot_table(
        index='group', columns='resolution',
        values=['ARI', 'NMI', 'homogeneity', 'completeness'],
    )
    wide.to_csv(outdir / 'metrics_wide.tsv', sep='\t')
    print(f'Saved wide-format pivot: {outdir/"metrics_wide.tsv"}')

    # Cluster N absorption across resolutions
    c_df = None
    if target_breakdown:
        c_df = (pd.DataFrame(target_breakdown).fillna(0).astype(int)
                  .sort_index())
        c_df.index.name = 'leiden_cluster'
        c_df.columns.name = 'resolution'
        out_path = outdir / f'cluster{args.target_cluster}_absorption.tsv'
        c_df.to_csv(out_path, sep='\t')
        print(f'Saved cluster {args.target_cluster} absorption: {out_path}')
        print(f'\nCluster {args.target_cluster} absorption across resolutions:')
        print(c_df.to_string())

    # Optional figures
    if args.plots:
        print(f'\nRendering figures ({", ".join(args.plot_format)}):')
        plot_metrics_sweep(master, outdir, args.plot_format)
        if c_df is not None:
            plot_target_absorption(c_df, outdir, args.plot_format,
                                   args.target_cluster)

    print('\nDone.')


if __name__ == '__main__':
    main()
