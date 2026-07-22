#!/usr/bin/env python3
"""
paga_trajectory.py — Fig 4A: PAGA connectivity graph over the 18 clusters.

Cleaned, parameterised port of Z/paga_iteration.ipynb (the original iteration
sweep over n_neighbors x threshold). The notebook's final left-over settings were
n_neighbors=8, n_pcs=20, threshold=0.09; the manuscript Fig 4A caption records
edge threshold 0.3 — CONFIRM which is canonical and pass via --threshold.

PAGA nodes are coloured with the SAME 18-shade palette (data/cluster_colors_18.json)
and ordered by data/cluster_order.json (adeno -> NE phenotypic axis), so cluster
colours match every other figure. Cluster 18 (inferCNV reference) is dropped.

Reads an object that already carries X_pca / X_umap / seurat_clusters (the original
Seurat-exported h5ad). Recomputes neighbors before PAGA (as the notebook did).

Example:
  python paga_trajectory.py --h5ad finaldong-log.h5ad \\
      --cluster-key seurat_clusters --n-neighbors 8 --n-pcs 20 --threshold 0.3 \\
      --cluster-order data/cluster_order.json --colors data/cluster_colors_18.json \\
      --outdir figures --stem fig4A_paga
"""
import argparse, json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--n-neighbors', type=int, default=8, help='number of neightbors for the drawn graph = 8')
    ap.add_argument('--n-pcs', type=int, default=20, help='pcs for the drawn graph = 20')
    ap.add_argument('--threshold', type=float, default=0.09,
                    help='edge threshold for the drawn graph = 0.09')
    ap.add_argument('--recompute-pca', action='store_true',
                    help='re-run sc.tl.pca(svd_solver="arpack") as the notebook did')
    ap.add_argument('--flip', choices=['none', 'v', 'h', 'both'], default='none',
                    help='reflect the PAGA layout: v=vertical (fixes an upside-down graph), '
                         'h=horizontal, both=180deg. FR layout orientation is otherwise arbitrary')
    ap.add_argument('--random-state', type=int, default=0,
                    help='layout seed; FR orientation is otherwise non-deterministic run-to-run')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--colors', default='data/cluster_colors_18.json',
                    help='flat {cluster: hex} 18-shade palette')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='cluster ids to drop (default 18 = inferCNV reference)')
    ap.add_argument('--compare', action='store_true',
                    help='also write the paga_compare panel (UMAP beside graph)')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig4A_paga')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa')
    ap.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm')
    args = ap.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        w = (args.figwidth_mm / 25.4) * w_mult if args.figwidth_mm else default_w_in
        if args.figheight_mm:
            h = (args.figheight_mm / 25.4) * h_mult
        elif args.figwidth_mm:
            h = w * (default_h_in / default_w_in)  # width-only: preserve authored aspect (avoid skew)
        else:
            h = default_h_in
        return (w, h)

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['font.size'] = 10       # cluster legend/colorbar text
    matplotlib.rcParams['axes.titlesize'] = 13
    matplotlib.rcParams['legend.fontsize'] = 8
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc

    sc.settings.verbosity = 1
    sc.settings.figdir = args.outdir
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            print(f'excluding clusters {args.exclude_clusters}: dropping {int((~keep).sum())} cells')
            adata = adata[keep].copy()

    # clusters as an ordered categorical following the canonical phenotypic axis,
    # so PAGA node colours line up with cluster_colors_18.json
    clusters = adata.obs[args.cluster_key].astype(str)
    present = sorted(clusters.unique(), key=lambda s: int(s) if s.isdigit() else 1e9)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        cl_order = [c for c in co if c in present] + [c for c in present if c not in co]
    else:
        cl_order = present
    adata.obs['clusters'] = clusters.astype('category').cat.reorder_categories(cl_order)

    if args.colors:
        palette = {str(k): v for k, v in json.load(open(args.colors)).items()}
        adata.uns['clusters_colors'] = [palette.get(c, '#cccccc') for c in cl_order]
        missing = [c for c in cl_order if c not in palette]
        if missing:
            print(f'[warn] no colour for clusters {missing} -> grey')

    if args.recompute_pca or 'X_pca' not in adata.obsm:
        print('computing PCA (arpack)')
        sc.tl.pca(adata, svd_solver='arpack')

    print(f'neighbors: n_neighbors={args.n_neighbors}, n_pcs={args.n_pcs}')
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.paga(adata, groups='clusters')

    # Compute the layout ONCE (deterministic seed), then reflect the node positions so the
    # orientation is fixed and identical across every output / paga_compare. FR layout
    # otherwise picks an arbitrary rotation/flip each run -> the "upside-down" graph.
    tmpfig, tmpax = plt.subplots()
    sc.pl.paga(adata, threshold=args.threshold, color='clusters',
               random_state=args.random_state, ax=tmpax, show=False)
    plt.close(tmpfig)
    pos = np.asarray(adata.uns['paga']['pos'], dtype=float)
    if args.flip in ('v', 'both'):
        pos[:, 1] = -pos[:, 1]
    if args.flip in ('h', 'both'):
        pos[:, 0] = -pos[:, 0]
    if args.flip != 'none':
        print(f'flipping layout: {args.flip}')

    # --- PAGA graph (Fig 4A) ---
    for fmt in args.plot_format:
        fig, ax = plt.subplots(figsize=figsize(6, 6))
        sc.pl.paga(adata, threshold=args.threshold, color='clusters', pos=pos,
                   node_size_scale=0.6, fontsize=8, fontoutline=1,#fontsize=11, fontoutline=1,
                   frameon=False, ax=ax, show=False)
        ax.set_axis_off()   # belt-and-suspenders: also drop any ticks scanpy leaves behind
        fig.savefig(Path(args.outdir) / f'{args.stem}_thr{args.threshold}.{fmt}',
                    bbox_inches='tight', dpi=300)
        plt.close(fig)
    print(f'[ok] wrote {args.stem}_thr{args.threshold}.{{{",".join(args.plot_format)}}}')

    # --- optional paga_compare (graph beside UMAP) ---
    if args.compare:
        if 'X_umap' not in adata.obsm:
            print('[warn] no X_umap; skipping paga_compare')
        else:
            for fmt in args.plot_format:
                axs = sc.pl.paga_compare(adata, threshold=args.threshold, pos=pos,
                                         node_size_scale=0.6, fontsize=11, fontoutline=1,
                                         show=False)
                fig = axs[0].figure if isinstance(axs, (list, tuple)) else axs.figure
                fig.savefig(Path(args.outdir) / f'{args.stem}_compare_thr{args.threshold}.{fmt}',
                            bbox_inches='tight', dpi=300)
                plt.close(fig)
            print(f'[ok] wrote {args.stem}_compare_thr{args.threshold}.*')


if __name__ == '__main__':
    main()
