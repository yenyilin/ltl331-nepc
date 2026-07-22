#!/usr/bin/env python3
"""
cluster_origin_correspondence.py

Unsupervised hierarchical clustering of the origin-vs-Harmony row-fraction matrix
produced by plot_origin_vs_harmony.py (rows = dataset_originalCluster, cols =
Harmony joint cluster, values = fraction of that original cluster's cells in each
Harmony cluster).

Clustering the ROWS groups original clusters that distribute into the same Harmony
clusters -> cross-dataset correspondence. Original clusters from different
datasets falling in one group are LTL331 <-> clinical (Gao/Li) matches, recovered
without any supervision.

Outputs in OUTDIR:
  <prefix>_clustermap.{pdf,png}    clustermap with dataset colour bar (fonttype 42)
  <prefix>_groups.tsv              every original cluster -> group assignment
  <prefix>_crossdataset.tsv        only the groups spanning >1 dataset

Usage:
  python cluster_origin_correspondence.py \\
      --rowfrac origin_map_055/origin_vs_harmony_rowfrac.tsv \\
      --outdir origin_corr_055 \\
      --metric cosine --method average \\
      --criterion distance --t 0.5 \\
      --plot-format pdf png
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from _utils import setup_style, numkey


def dataset_of(label, ltl331_label):
    ds = label.rsplit('_', 1)[0]
    return 'LTL331' if ds in (ltl331_label, 'other', 'LTL331') else ds


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--rowfrac', required=True,
                   help='origin_vs_harmony_rowfrac.tsv')
    p.add_argument('--outdir', default='origin_corr')
    p.add_argument('--metric', default='cosine',
                   choices=['cosine', 'correlation', 'euclidean', 'jensenshannon'])
    p.add_argument('--method', default='average',
                   help='scipy linkage method (average, complete, ward, ...)')
    p.add_argument('--criterion', default='distance',
                   choices=['distance', 'maxclust'])
    p.add_argument('--t', type=float, default=0.5,
                   help='distance cut (criterion=distance) or n_clusters '
                        '(criterion=maxclust)')
    p.add_argument('--ltl331-label', default='LTL331',
                   help='dataset prefix for LTL331 rows (also maps "other")')
    p.add_argument('--prefix', default='origin_correspondence')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    frac = pd.read_csv(args.rowfrac, sep='\t', index_col=0).fillna(0)
    frac = frac[sorted(frac.columns, key=numkey)]
    print(f'Loaded {frac.shape[0]} original clusters x {frac.shape[1]} '
          f'Harmony clusters')

    ds = pd.Series([dataset_of(l, args.ltl331_label) for l in frac.index],
                   index=frac.index, name='dataset')

    # linkage (ward needs euclidean)
    method = args.method
    metric = args.metric
    if method == 'ward' and metric != 'euclidean':
        print('[warn] ward requires euclidean; switching metric to euclidean')
        metric = 'euclidean'
    Z = linkage(pdist(frac.values, metric=metric), method=method)

    if args.criterion == 'maxclust':
        groups = fcluster(Z, t=int(args.t), criterion='maxclust')
    else:
        groups = fcluster(Z, t=args.t, criterion='distance')

    tab = pd.DataFrame({'original': frac.index, 'dataset': ds.values,
                        'group': groups}).sort_values(
                            ['group', 'dataset', 'original'])
    tab.to_csv(outdir / f'{args.prefix}_groups.tsv', sep='\t', index=False)

    print(f'\n=== groups (metric={metric}, method={method}, '
          f'{args.criterion}={args.t}) ===')
    print(f'total groups: {tab["group"].nunique()}')
    print('\ncross-dataset groups (span >1 dataset):')
    cross_rows = []
    for gid, sub in tab.groupby('group'):
        if sub['dataset'].nunique() > 1:
            datasets = sorted(sub['dataset'].unique())
            print(f'  group {gid} [{"+".join(datasets)}]: '
                  f'{sub["original"].tolist()}')
            for m, d in zip(sub['original'], sub['dataset']):
                cross_rows.append({'group': gid, 'datasets': '+'.join(datasets),
                                   'original': m, 'dataset': d})
    if cross_rows:
        pd.DataFrame(cross_rows).to_csv(
            outdir / f'{args.prefix}_crossdataset.tsv', sep='\t', index=False)
    else:
        print('  (none at this threshold — try a larger --t or '
              'criterion=maxclust)')

    # clustermap
    setup_style(font_size=7)
    import seaborn as sns
    from matplotlib.patches import Patch

    palette = {'LTL331': 'tab:blue', 'Gao': 'tab:orange', 'Li': 'tab:green'}
    row_colors = ds.map(palette)

    g = sns.clustermap(
        frac, row_linkage=Z, col_cluster=False,
        row_colors=row_colors, cmap='viridis', vmin=0, vmax=1,
        figsize=(frac.shape[1] * 0.4 + 4, frac.shape[0] * 0.22 + 4),
        yticklabels=True, xticklabels=True,
        cbar_kws={'label': 'frac of original cluster'},
    )
    g.ax_heatmap.set_xlabel('Harmony joint cluster')
    g.ax_heatmap.set_ylabel('original cluster')
    handles = [Patch(facecolor=c, label=d) for d, c in palette.items()]
    g.ax_heatmap.legend(handles=handles, title='dataset',
                        bbox_to_anchor=(1.28, 1.0), loc='upper left',
                        frameon=False, fontsize=6)

    for fmt in args.plot_format:
        path = outdir / f'{args.prefix}_clustermap.{fmt}'
        g.savefig(path, format=fmt)
        print(f'wrote {path}')

    print(f'\nSaved group tables to {outdir}/')
    print('Done.')


if __name__ == '__main__':
    main()
