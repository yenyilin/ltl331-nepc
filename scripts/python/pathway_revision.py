#!/usr/bin/env python3
"""
pathway_revision.py — disambiguate the KRAS_UP / KRAS_DN co-enrichment in
NE clusters AND analyse per-timepoint oncogenic-pathway dynamics in cluster 4
(survival state) and cluster 12 (MYC).

Three analyses, all driven by per-cell gene-set scoring (scanpy.tl.score_genes,
identical to sc.score_genes — matched-background; deterministic with --seed):

  (1) Per-cluster mean scores for every gene set
      -> per_cluster_pathway_scores.tsv + heatmap
  (2) KRAS pair disambiguation (--kras-pair UP_name DN_name):
      - Jaccard overlap of UP vs DN gene sets
      - per-cluster UP, DN and directional (UP - DN) scores
      - top driver genes per set per cluster (which genes are doing the "lift")
      -> kras_directional.tsv + plot, kras_top_drivers.tsv
  (3) Per-(target cluster x timepoint) pathway scores (--target-clusters)
      -> cluster<N>_timepoint_pathway_scores.tsv + line plots

Input gene sets file (TSV):
    pathway_name<TAB>gene1,gene2,...,geneN

Get the Hallmark KRAS_SIGNALING_UP/DN gene lists from MSigDB and place them in
the TSV (you can include any Hallmark / C2:CGP sets you want scored).

Deps: scanpy/anndata, numpy, pandas, matplotlib, seaborn (optional).

Example:
  python pathway_revision.py \\
      --h5ad ltl331.h5ad --cluster-key clusters \\
      --gene-sets hallmark_sets.tsv \\
      --kras-pair HALLMARK_KRAS_SIGNALING_UP HALLMARK_KRAS_SIGNALING_DN \\
      --target-clusters 4 12 --timepoint-key sample \\
      --timepoint-order preCx 2wk 4wk 8wk 12wk 16wk 20wk 22wk \\
      --outdir pathway_revision --plot-format pdf png
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _utils import setup_style, savefig, timepoint_key


def load_gene_sets(path):
    sets = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        name, genes = line.split('\t', 1)
        sets[name.strip()] = [g.strip() for g in genes.replace(',', ' ').split()
                              if g.strip()]
    return sets


def score_sets(adata, sets, seed, use_raw=False):
    """scanpy.tl.score_genes for each set; returns DataFrame of per-cell scores."""
    import scanpy as sc
    var_names = adata.raw.var_names if use_raw else adata.var_names
    scores = {}
    for name, genes in sets.items():
        present = [g for g in genes if g in var_names]
        if len(present) < 3:
            print(f'[score] {name}: only {len(present)} genes present, skipping')
            continue
        key = f'_score_{name}'
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=seed, use_raw=use_raw)
        scores[name] = adata.obs[key].values
        adata.obs.drop(columns=[key], inplace=True)
        print(f'[score] {name}: {len(present)}/{len(genes)} genes used')
    return pd.DataFrame(scores, index=adata.obs_names)


def per_cluster_means(scores_df, clusters):
    s = scores_df.copy()
    s['cluster'] = clusters.values
    return s.groupby('cluster').mean()


def jaccard(a, b):
    A, B = set(a), set(b)
    u = A | B
    return len(A & B) / len(u) if u else 0.0


def top_driver_genes(adata, set_genes, cluster_mask, k=10, use_raw=False):
    """Top genes in a set ranked by mean expression in a given cluster."""
    import scipy.sparse as sp
    src = adata.raw if use_raw else adata
    genes = [g for g in set_genes if g in src.var_names]
    if not genes:
        return []
    idx = [src.var_names.get_loc(g) for g in genes]
    X = src.X[:, idx][cluster_mask]
    means = np.asarray(X.mean(axis=0)).ravel() if sp.issparse(X) \
            else np.asarray(X.mean(axis=0)).ravel()
    order = np.argsort(means)[::-1][:k]
    return list(zip(np.array(genes)[order], np.round(means[order], 3)))


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--gene-sets', required=True,
                   help='TSV: pathway_name<TAB>comma-separated-genes')
    p.add_argument('--cluster-key', default='clusters')
    p.add_argument('--kras-pair', nargs=2, default=None,
                   metavar=('UP_NAME', 'DN_NAME'),
                   help='names of the two paired gene sets for directional KRAS')
    p.add_argument('--target-clusters', nargs='+', default=None,
                   help='clusters to score across timepoints (e.g. 4 12)')
    p.add_argument('--timepoint-key', default=None)
    p.add_argument('--timepoint-order', nargs='+', default=None)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--use-raw', action='store_true',
                   help='score from adata.raw (lognorm, full gene set) instead of '
                        '.X; use when .X holds Seurat scaled/z-scored data')
    p.add_argument('--outdir', default='pathway_revision')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    try:
        import anndata as ad
    except Exception as e:
        raise SystemExit(f'anndata required ({e}); activate scanpy env')
    adata = ad.read_h5ad(args.h5ad)
    clusters = adata.obs[args.cluster_key].astype(str)
    sets = load_gene_sets(args.gene_sets)
    print(f'Loaded {adata.n_obs} cells; {len(sets)} gene sets')

    # --- (1) per-cell scores + per-cluster means ---
    scores = score_sets(adata, sets, args.seed, use_raw=args.use_raw)
    if scores.empty:
        raise SystemExit('no scorable gene sets')
    per_cl = per_cluster_means(scores, clusters)
    per_cl.to_csv(outdir / 'per_cluster_pathway_scores.tsv', sep='\t')
    print('\n[1] per-cluster pathway score summary:')
    print(per_cl.round(3).to_string())

    plt = setup_style()
    try:
        import seaborn as sns
        order = sorted(per_cl.index, key=lambda c: int(c) if c.isdigit() else 1e9)
        z = per_cl.loc[order]
        z = (z - z.mean()) / z.std().replace(0, np.nan)
        fig = plt.figure(figsize=(max(6, z.shape[1] * 0.5 + 2),
                                  max(4, z.shape[0] * 0.25 + 2)))
        sns.heatmap(z.fillna(0), cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'z-score'})
        plt.xticks(rotation=45, ha='right', fontsize=6)
        plt.title('Per-cluster pathway scores (z across clusters)')
        plt.tight_layout()
        savefig(fig, outdir, 'per_cluster_pathway_heatmap', args.plot_format)
        plt.close(fig)
    except Exception as e:
        print(f'[1] heatmap skipped: {e}')

    # --- (2) KRAS directional ---
    if args.kras_pair:
        up, dn = args.kras_pair
        if up not in sets or dn not in sets:
            print(f'[2] KRAS pair: missing set(s) in gene-sets file: {up}/{dn}')
        elif up not in scores.columns or dn not in scores.columns:
            print(f'[2] KRAS pair: not enough genes for one of {up}/{dn}')
        else:
            j = jaccard(sets[up], sets[dn])
            print(f'\n[2] {up} vs {dn} gene-set Jaccard overlap: {j:.3f}')
            shared = sorted(set(sets[up]) & set(sets[dn]))
            (pd.DataFrame({'shared_genes': shared})
               .to_csv(outdir / 'kras_set_overlap.tsv', sep='\t', index=False))

            tab = per_cl[[up, dn]].copy()
            tab['directional_UP_minus_DN'] = tab[up] - tab[dn]
            tab = tab.sort_values('directional_UP_minus_DN', ascending=False)
            tab.to_csv(outdir / 'kras_directional_per_cluster.tsv', sep='\t')
            print(tab.round(3).to_string())

            # plot: directional bar + side panel of UP and DN means
            cl_order = list(tab.index)
            fig, axes = plt.subplots(1, 2, figsize=(max(10, len(cl_order)*0.45+4), 3))
            axes[0].bar(range(len(cl_order)), tab['directional_UP_minus_DN'].values,
                        color=['#1f77b4' if v < 0 else '#d62728'
                               for v in tab['directional_UP_minus_DN']])
            axes[0].axhline(0, color='k', lw=0.4)
            axes[0].set_xticks(range(len(cl_order)))
            axes[0].set_xticklabels(cl_order, fontsize=7)
            axes[0].set_ylabel('UP − DN score')
            axes[0].set_title('Directional KRAS (UP − DN) per cluster')
            axes[1].plot(range(len(cl_order)), tab[up].values, 'o-', label=up,
                         color='#d62728')
            axes[1].plot(range(len(cl_order)), tab[dn].values, 'o-', label=dn,
                         color='#1f77b4')
            axes[1].set_xticks(range(len(cl_order)))
            axes[1].set_xticklabels(cl_order, fontsize=7)
            axes[1].set_ylabel('mean score'); axes[1].legend(fontsize=6, frameon=False)
            axes[1].set_title('Raw UP / DN means per cluster')
            savefig(fig, outdir, 'kras_directional', args.plot_format)
            plt.close(fig)

            # top driver genes per cluster (which genes drive each set's signal)
            ne_clusters = args.target_clusters or sorted(per_cl.index,
                                                        key=lambda c: int(c) if c.isdigit() else 1e9)
            rows = []
            for c in ne_clusters:
                mask = (clusters.values == str(c))
                if mask.sum() == 0:
                    continue
                for set_name, set_label in ((up, 'UP'), (dn, 'DN')):
                    drivers = top_driver_genes(adata, sets[set_name], mask, k=10, use_raw=args.use_raw)
                    for rank, (g, mn) in enumerate(drivers, 1):
                        rows.append({'cluster': c, 'set': set_label, 'rank': rank,
                                     'gene': g, 'mean_in_cluster': mn})
            pd.DataFrame(rows).to_csv(outdir / 'kras_top_drivers.tsv',
                                       sep='\t', index=False)
            print(f'[2] wrote top driver genes for clusters: {ne_clusters}')

    # --- (3) per-timepoint pathway scores in target clusters ---
    if args.target_clusters and args.timepoint_key:
        tp = adata.obs[args.timepoint_key].astype(str)
        if args.timepoint_order:
            order = [str(t) for t in args.timepoint_order]
        else:
            order = sorted(tp.unique(), key=timepoint_key)

        scores_df = scores.copy()
        scores_df['cluster'] = clusters.values
        scores_df['timepoint'] = tp.values
        for c in args.target_clusters:
            sub = scores_df[scores_df['cluster'] == str(c)]
            if sub.empty:
                print(f'[3] cluster {c}: no cells')
                continue
            per_t = sub.drop(columns=['cluster']).groupby('timepoint').mean()
            # Do NOT forward-fill: a timepoint where the cluster has no cells must
            # read as a gap (n=0), not a copy of the previous timepoint's scores.
            n_per_t = sub.groupby('timepoint').size().reindex(order).fillna(0).astype(int)
            per_t = per_t.reindex(order)
            per_t.insert(0, 'n_cells', n_per_t.values)
            per_t.to_csv(outdir / f'cluster{c}_timepoint_pathway_scores.tsv',
                         sep='\t')
            print(f'\n[3] cluster {c} per-timepoint pathway scores:')
            print(per_t.round(3).to_string())

            fig, ax = plt.subplots(figsize=(max(6, len(order) * 0.7 + 2), 3.5))
            for col in per_t.columns:
                if col == 'n_cells':
                    continue
                ax.plot(range(len(order)), per_t[col].values, marker='o',
                        linewidth=1.2, label=col)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, fontsize=8)
            ax.set_ylabel('pathway score (per-cell mean within cluster)')
            ax.set_title(f'Cluster {c} pathway dynamics')
            ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=6,
                      frameon=False)
            savefig(fig, outdir, f'cluster{c}_pathway_dynamics', args.plot_format)
            plt.close(fig)

    print('\nDone.')


if __name__ == '__main__':
    main()
