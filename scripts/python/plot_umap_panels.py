#!/usr/bin/env python3
"""
plot_umap_panels.py — Fig 2A (UMAP coloured by cluster) and Fig 2B (feature UMAPs, e.g.
AR / NCAM1) from the h5ad. Self-contained matplotlib (full control, small rasterized PDFs).

  --which clusters   F2A: UMAP coloured by seurat_clusters using the 18-shade phenotype
                     palette (data/cluster_colors_18.json) so clusters group visually by
                     phenotype; cluster IDs labelled at centroids + a legend (id -> phenotype).
  --which features   F2B: one UMAP per --features gene, coloured by expression (from .raw
                     lognorm — .X is scaled), high-expressing cells drawn on top.
  --which both       both (default).

CRITICAL: the LTL331 base object's .X is SCALED; expression for F2B is read from .raw.

Examples:
  python plot_umap_panels.py --h5ad ltl331_base.h5ad --which clusters \\
      --cluster-colors data/cluster_colors_18.json --labels data/cluster_phenotype_labels.json \\
      --cluster-order data/cluster_order.json --outdir figures --stem fig2

  python plot_umap_panels.py --h5ad ltl331_base.h5ad --which features \\
      --features AR NCAM1 --outdir figures --stem fig2
"""
import argparse, json
from pathlib import Path

import numpy as np
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--which', choices=['clusters', 'features', 'both', 'categorical'],
                    default='both')
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--umap-key', default='X_umap')
    # generic categorical mode (Fig 6A joint cohort/lineage; Fig 2 timepoint panel)
    ap.add_argument('--color-by', default=None,
                    help='obs categorical column to colour by (--which categorical)')
    ap.add_argument('--facet-by', default=None,
                    help='obs column to split into subpanels (--which categorical)')
    ap.add_argument('--color-map', default=None,
                    help='optional json {value: hex} palette for --color-by')
    ap.add_argument('--na-color', default='#dddddd',
                    help='grey context colour for faceted background')
    ap.add_argument('--cluster-colors', default='data/cluster_colors_18.json')
    ap.add_argument('--labels', default='data/cluster_phenotype_labels.json',
                    help='cluster->label json for the F2A legend (optional)')
    ap.add_argument('--cluster-order', default='data/cluster_order.json',
                    help='legend order (optional)')
    ap.add_argument('--features', nargs='+', default=['AR', 'NCAM1'])
    ap.add_argument('--feature-ncols', type=int, default=None,
                    help='columns for the F2B feature grid (default: all in one row). '
                         'e.g. 1 stacks features vertically (2 features -> 2 rows).')
    ap.add_argument('--use-raw', dest='use_raw', action='store_true', default=True)
    ap.add_argument('--no-use-raw', dest='use_raw', action='store_false')
    ap.add_argument('--cmap', default='viridis', help='feature-UMAP colormap (F2B)')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'],
                    help='clusters to drop (default 18 = inferCNV reference)')
    ap.add_argument('--point-size', type=float, default=3.0)
    ap.add_argument('--label-centroids', dest='label_centroids', action='store_true', default=True)
    ap.add_argument('--no-label-centroids', dest='label_centroids', action='store_false')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig2')
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author the panel at this FINAL width in mm (the assemble_figure --qa '
                         'onpage_mm target) so the assembler scales it ~1:1 and fonts stay readable; '
                         'the saved PDF is wider than this by any outside legend — tune with --qa. '
                         'This is the TOTAL panel width (facets divide it); height follows '
                         'the authored aspect.')
    ap.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        """final size in inches; --figwidth-mm/--figheight-mm override the inch defaults.
        w_mult/h_mult scale the mm width/height for multi-facet grids (per-facet mm)."""
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh:
            return ((fw / 25.4) * w_mult, (fh / 25.4) * h_mult)
        if fw:   # width only -> preserve authored aspect ratio (avoids skew + empty margins)
            w = (fw / 25.4) * w_mult
            return (w, default_h_in * w / default_w_in)
        if fh:   # height only -> preserve authored aspect ratio
            h = (fh / 25.4) * h_mult
            return (default_w_in * h / default_h_in, h)
        return (default_w_in, default_h_in)

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    if args.exclude_clusters and args.cluster_key in adata.obs:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            print(f'excluding clusters {args.exclude_clusters}: dropping {int((~keep).sum())} cells')
            adata = adata[keep.values].copy()
    if args.umap_key not in adata.obsm:
        raise SystemExit(f'{args.umap_key} not in adata.obsm ({list(adata.obsm)})')
    um = np.asarray(adata.obsm[args.umap_key])[:, :2]
    n = um.shape[0]
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    def _csort(vals):
        return sorted(set(vals), key=lambda s: (0, int(s)) if str(s).isdigit() else (1, str(s)))

    # ---- generic categorical UMAP (Fig 6A joint cohort/lineage; Fig 2 timepoint) ----
    if args.which == 'categorical':
        if not args.color_by or args.color_by not in adata.obs:
            raise SystemExit(f'--which categorical needs --color-by <obs col>; got '
                             f'{args.color_by!r} (obs: {list(adata.obs)[:15]}...)')
        cval = adata.obs[args.color_by].astype(str).values
        cats = _csort(cval)
        if args.color_map and Path(args.color_map).exists():
            pal = {str(k): v for k, v in json.load(open(args.color_map)).items()}
        else:
            base = plt.get_cmap('tab10' if len(cats) <= 10 else 'tab20')
            pal = {c: base(i % base.N) for i, c in enumerate(cats)}
        legend = [Line2D([0], [0], marker='o', lw=0, markersize=6, markeredgecolor='none',
                         markerfacecolor=pal.get(c, '#999'), label=c) for c in cats]

        if args.facet_by:
            if args.facet_by not in adata.obs:
                raise SystemExit(f'--facet-by {args.facet_by!r} not in obs')
            fval = adata.obs[args.facet_by].astype(str).values
            facets = _csort(fval)
            fig, axes = plt.subplots(1, len(facets),
                                     figsize=figsize(4.2 * len(facets), 4.2),
                                     squeeze=False)
            for j, fc in enumerate(facets):
                ax = axes[0][j]
                ax.scatter(um[:, 0], um[:, 1], s=args.point_size * 0.5, c=args.na_color,
                           lw=0, rasterized=True)                      # grey context = all cells
                m = fval == fc
                sh = rng.permutation(int(m.sum()))
                xm, ym, cm = um[m, 0], um[m, 1], cval[m]
                ax.scatter(xm[sh], ym[sh], s=args.point_size,
                           c=[pal.get(cm[i], '#999') for i in sh], lw=0, rasterized=True)
                ax.set_title(f'{fc} (n={int(m.sum()):,})', fontsize=9)
                ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP1', fontsize=8)
                if j == 0:
                    ax.set_ylabel('UMAP2', fontsize=8)
                for s in ax.spines.values():
                    s.set_visible(False)
            axes[0][-1].legend(handles=legend, loc='center left', bbox_to_anchor=(1.005, 0.5),
                               fontsize=7, frameon=False, title=args.color_by)
            suffix = f'_{args.color_by}_by_{args.facet_by}'
        else:
            fig, ax = plt.subplots(figsize=figsize(6.5, 5.5))
            sh = rng.permutation(n)
            ax.scatter(um[sh, 0], um[sh, 1], s=args.point_size,
                       c=[pal.get(cval[i], '#999') for i in sh], lw=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.legend(handles=legend, loc='center left', bbox_to_anchor=(1.005, 0.5),
                      fontsize=7, frameon=False, title=args.color_by)
            suffix = f'_{args.color_by}'
        for fmt in args.plot_format:
            fig.savefig(outdir / f'{args.stem}{suffix}.{fmt}', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'[ok] wrote {args.stem}{suffix}.{{{",".join(args.plot_format)}}} '
              f'(colour={args.color_by}'
              + (f', facet={args.facet_by}' if args.facet_by else '') + ')')
        return

    # ---- F2A: UMAP by cluster ----
    if args.which in ('clusters', 'both'):
        cl = adata.obs[args.cluster_key].astype(str).values
        cats = sorted(set(cl), key=lambda s: int(s) if s.isdigit() else 1e9)
        colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
            if Path(args.cluster_colors).exists() else {}
        if not colors:
            tab = plt.get_cmap('tab20')
            colors = {c: tab(i % 20) for i, c in enumerate(cats)}
        labels = {str(k): v for k, v in json.load(open(args.labels)).items()} \
            if (args.labels and Path(args.labels).exists()) else {}
        if args.cluster_order and Path(args.cluster_order).exists():
            co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
            leg_order = [c for c in co if c in cats] + [c for c in cats if c not in co]
        else:
            leg_order = cats

        fig, ax = plt.subplots(figsize=figsize(6.5, 5.5))
        sh = rng.permutation(n)               # shuffle so no cluster sits permanently on top
        ax.scatter(um[sh, 0], um[sh, 1], s=args.point_size,
                   c=[colors.get(cl[i], '#cccccc') for i in sh],
                   lw=0, rasterized=True)
        if args.label_centroids:
            for c in cats:
                m = cl == c
                cx, cy = np.median(um[m, 0]), np.median(um[m, 1])
                ax.text(cx, cy, c, fontsize=8, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
        for s in ax.spines.values():
            s.set_visible(False)
        handles = [Line2D([0], [0], marker='o', lw=0, markersize=6,
                          markerfacecolor=colors.get(c, '#cccccc'), markeredgecolor='none',
                          label=labels.get(c, c)) for c in leg_order]
        ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.005, 0.5),
                  fontsize=7, frameon=False, handletextpad=0.3, labelspacing=0.35)
        ax.set_title(f'{n:,} transcriptomes — clusters', fontsize=9)
        for fmt in args.plot_format:
            fig.savefig(outdir / f'{args.stem}_clusters.{fmt}', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'[ok] wrote {args.stem}_clusters.{{{",".join(args.plot_format)}}}')

    # ---- F2B: feature UMAPs ----
    if args.which in ('features', 'both'):
        src = adata.raw if (args.use_raw and adata.raw is not None) else adata
        var = list(map(str, src.var_names)); idx = {g: i for i, g in enumerate(var)}
        feats = [g for g in args.features if g in idx]
        miss = [g for g in args.features if g not in idx]
        if miss:
            print(f'[warn] features absent from {"raw" if src is adata.raw else ".X"}: {miss}')
        if not feats:
            raise SystemExit('no requested features present in the object')
        k = len(feats)
        ncols = min(args.feature_ncols or k, k)   # default: single row
        nrows = int(np.ceil(k / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=figsize(4.6 * ncols, 4.4 * nrows),
                                 squeeze=False)
        for j, g in enumerate(feats):
            ax = axes[j // ncols][j % ncols]
            col = src.X[:, idx[g]]
            e = (np.asarray(col.todense()).ravel() if sp.issparse(col)
                 else np.asarray(col).ravel())
            o = np.argsort(e)                 # low first, high drawn on top
            sca = ax.scatter(um[o, 0], um[o, 1], s=args.point_size, c=e[o],
                             cmap=args.cmap, lw=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(g, fontsize=9, fontstyle='italic')
            ax.set_xlabel('UMAP1', fontsize=6)
            if j % ncols == 0:
                ax.set_ylabel('UMAP2', fontsize=6)
            for s in ax.spines.values():
                s.set_visible(False)
            cb = fig.colorbar(sca, ax=ax, shrink=0.6, aspect=12, pad=0.02)
            cb.set_label('lognorm expr.', fontsize=7); cb.ax.tick_params(labelsize=6)
        for e in range(k, nrows * ncols):     # hide any unused cells
            axes[e // ncols][e % ncols].axis('off')
        for fmt in args.plot_format:
            fig.savefig(outdir / f'{args.stem}_features.{fmt}', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'[ok] wrote {args.stem}_features.{{{",".join(args.plot_format)}}}  ({feats})')


if __name__ == '__main__':
    main()
