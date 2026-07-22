#!/usr/bin/env python3
"""
velocity_dynamical.py — Fig 4B (RNA-velocity stream embedding) and Fig 4C
(velocity pseudotime), rendered with the same 18-shade cluster palette as every
other figure.

Velocity is estimated in scVelo DYNAMICAL mode (recover_dynamics -> velocity ->
velocity_graph), where we choose full per-gene kinetics for a non-equilibrium transition).
The velocity computation here mirrors the one in cellrank_bifurcation.py so the
embedding and the fate analysis use identical velocities. If the object already
carries a `velocity` layer (and velocity_graph), plotting is done directly;
pass --recompute-velocity to (re)estimate from spliced/unspliced layers.

  --which stream      4B: velocity_embedding_stream on the UMAP, coloured by cluster
  --which pseudotime  4C: velocity_pseudotime on the UMAP (viridis)
  --which both        both (default)

NOTE on direction: velocity DIRECTION is unreliable globally here (latent_time vs
experimental time rho=-0.47), so the trajectory result is anchored on
velocity_pseudotime (rho=+0.70), NOT on the velocity field. 4B is shown as the velocity
field over the embedding; the directional claims live in Fig 4D-G.

Example:
  python velocity_dynamical.py --h5ad ltl331_velocity.h5ad \\
      --cluster-key seurat_clusters --umap-key X_umap \\
      --cluster-order data/cluster_order.json --colors data/cluster_colors_18.json \\
      --which stream --outdir figures --stem fig4B_velocity
"""
import argparse, json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--umap-key', default='X_umap')
    ap.add_argument('--pseudotime-key', default='velocity_pseudotime',
                    help='obs key for Fig 4C; computed if absent')
    ap.add_argument('--recompute-velocity', action='store_true',
                    help='(re)estimate velocity from spliced/unspliced layers')
    ap.add_argument('--velocity-mode', choices=['dynamical', 'stochastic'], default='dynamical')
    ap.add_argument('--n-jobs', type=int, default=1)
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--colors', default='data/cluster_colors_18.json',
                    help='flat {cluster: hex} 18-shade palette')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='cluster ids to drop (default 18 = inferCNV reference)')
    ap.add_argument('--which', choices=['stream', 'pseudotime', 'both'], default='both')
    ap.add_argument('--legend-loc', default='right margin',
                    help="scVelo legend location ('on data' to label clusters in-plot)")
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig4B_velocity')
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
    matplotlib.rcParams['font.size'] = 10       # scVelo ticks/colorbar; cluster legend set per-call
    matplotlib.rcParams['axes.titlesize'] = 14
    matplotlib.rcParams['axes.labelsize'] = 12
    matplotlib.rcParams['legend.fontsize'] = 7
    import matplotlib.pyplot as plt
    import scanpy as sc
    import scvelo as scv

    scv.settings.verbosity = 1
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes')
    if args.cluster_key not in adata.obs:
        raise SystemExit(f'cluster key {args.cluster_key} not in obs')
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            print(f'excluding clusters {args.exclude_clusters}: dropping {int((~keep).sum())} cells')
            adata = adata[keep].copy()

    # cluster key -> ordered categorical along the canonical axis + matching colours
    clusters = adata.obs[args.cluster_key].astype(str)
    present = sorted(clusters.unique(), key=lambda s: int(s) if s.isdigit() else 1e9)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        cl_order = [c for c in co if c in present] + [c for c in present if c not in co]
    else:
        cl_order = present
    adata.obs[args.cluster_key] = clusters.astype('category').cat.reorder_categories(cl_order)
    if args.colors:
        palette = {str(k): v for k, v in json.load(open(args.colors)).items()}
        adata.uns[f'{args.cluster_key}_colors'] = [palette.get(c, '#cccccc') for c in cl_order]
        missing = [c for c in cl_order if c not in palette]
        if missing:
            print(f'[warn] no colour for clusters {missing} -> grey')

    # velocity: recompute (dynamical) or use what's stored
    if args.recompute_velocity:
        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
            raise SystemExit('spliced/unspliced layers required to recompute velocity')
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        if args.velocity_mode == 'dynamical':
            scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
            scv.tl.velocity(adata, mode='dynamical')
        else:
            scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
    else:
        if 'velocity' not in adata.layers:
            raise SystemExit('no velocity layer found; pass --recompute-velocity')
        if 'velocity_graph' not in adata.uns:
            print('velocity_graph not stored; computing it')
            scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)

    if args.umap_key not in adata.obsm:
        raise SystemExit(f'{args.umap_key} not in adata.obsm ({list(adata.obsm)})')
    basis = args.umap_key[2:] if args.umap_key.startswith('X_') else args.umap_key

    # ---- 4B: velocity stream embedding ----
    # NOTE: render a fresh figure per format. scVelo adds the legend/colorbar ("band")
    # as a lazy draw-time axis; reusing one figure across savefig calls with
    # bbox_inches='tight' crops the first format before that band is finalized, putting
    # it on the wrong side (PDF broke, PNG — saved second — looked fine).
    if args.which in ('stream', 'both'):
        for fmt in args.plot_format:
            fig, ax = plt.subplots(figsize=figsize(6, 6))
            scv.pl.velocity_embedding_stream(
                adata, basis=basis, color=args.cluster_key, ax=ax, show=False,
                legend_loc=args.legend_loc, title='RNA velocity (dynamical)',
                fontsize=12, legend_fontsize=7, size=20, alpha=0.6, linewidth=1.0, arrow_size=1.2)
            fig.canvas.draw()  # finalize draw-time band placement before tight crop
            fig.savefig(outdir / f'{args.stem}_stream.{fmt}', bbox_inches='tight', dpi=300)
            plt.close(fig)
        print(f'[ok] wrote {args.stem}_stream.{{{",".join(args.plot_format)}}}')

    # ---- 4C: velocity pseudotime ----
    if args.which in ('pseudotime', 'both'):
        if args.pseudotime_key not in adata.obs.columns:
            print(f'{args.pseudotime_key} not in obs; computing velocity_pseudotime')
            scv.tl.velocity_pseudotime(adata)
            ptk = 'velocity_pseudotime'
        else:
            ptk = args.pseudotime_key
        # Draw a NATIVE matplotlib colorbar rather than scVelo's. scVelo adds the
        # colorbar as a lazy draw-time artist whose gradient is anchored to pre-crop
        # figure coords; under bbox_inches='tight' the fill shifts left of its frame
        # (empty box on the right). A native fig.colorbar is a real figure axis that
        # transforms correctly on tight-crop — matching plot_fate_umap.py.
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        vmin = float(adata.obs[ptk].min())
        vmax = float(adata.obs[ptk].max())
        # fresh figure per format (same scVelo lazy-band reason as 4B above)
        for fmt in args.plot_format:
            fig, ax = plt.subplots(figsize=figsize(6, 6))
            scv.pl.scatter(adata, basis=basis, color=ptk, cmap='viridis', ax=ax, show=False,
                           fontsize=12, size=20, title='velocity pseudotime', colorbar=False)
            sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap='viridis')
            sm.set_array([])
            cb = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=12, pad=0.02)
            cb.set_label(ptk)
            fig.savefig(outdir / f'{args.stem}_pseudotime.{fmt}', bbox_inches='tight', dpi=300)
            plt.close(fig)
        print(f'[ok] wrote {args.stem}_pseudotime.{{{",".join(args.plot_format)}}}')


if __name__ == '__main__':
    main()
