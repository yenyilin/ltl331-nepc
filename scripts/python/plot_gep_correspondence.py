#!/usr/bin/env python3
"""
plot_gep_correspondence.py — show how the seven cNMF gene-expression programs (GEPs)
map onto the 18 seurat_clusters and the 8 timepoints.

Each GEP is scored per cell with sc.tl.score_genes over its top-N genes, averaged per
group (cluster / timepoint), then z-scored per GEP across groups so each row shows where
that program is relatively enriched. Rows (GEPs) are ordered by their peak cluster along
the canonical adenocarcinoma->neuroendocrine axis (data/cluster_order.json); the per-group
argmax is boxed so the "GEP4 -> AR-high PRAD, GEP6 -> ASCL1+ cl10" correspondence is read
directly off the diagonal.

Input GEP table (--gep): rows = rank-ordered genes, one column per GEP (data/gep_top100.tsv,
or the .xlsx — needs openpyxl). Genes are matched case-insensitively to adata.var_names.

Layer: score_genes needs per-cell expression. By default uses .raw if present, else .X;
pass --use-raw / --layer to force. If .X looks scaled (has sizeable negatives) and no raw
exists, scoring still ranks groups correctly (the per-GEP z-score absorbs the scale) but a
warning is printed.

Example:
  python plot_gep_correspondence.py --h5ad ltl331.h5ad --gep data/gep_top100.tsv \\
      --cluster-key seurat_clusters --timepoint-key timepoint --week-key week_num \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --outdir figures --stem fig_gep_correspondence
"""
import argparse, json, re
from pathlib import Path

import numpy as np
import pandas as pd

# cosmetic short labels (Table 1); override is not needed, purely for row annotation
GEP_LABELS = {
    'GEP1': 'NE ASCL1− (terminal)',
    'GEP2': 'AR-low luminal (castration-adapted)',
    'GEP3': 'EMT / stem-like (AR-low)',
    'GEP4': 'AR+ luminal (PRAD)',
    'GEP5': 'Proliferation (cell cycle)',
    'GEP6': 'NE ASCL1+',
    'GEP7': 'Ciliogenesis',
}


def load_geps(path, top):
    p = Path(path)
    if p.suffix.lower() in ('.xlsx', '.xls'):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep='\t')
    df = df[[c for c in df.columns if not str(c).startswith('Unnamed')]]
    geps = {}
    for c in df.columns:
        genes = [str(g).strip() for g in df[c].dropna().tolist() if str(g).strip()]
        geps[c] = genes[:top]
    return geps


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--gep', default='data/gep_top100.tsv')
    ap.add_argument('--top', type=int, default=100, help='top-N genes per GEP (default 100)')
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--timepoint-key', default='timepoint')
    ap.add_argument('--week-key', default='week_num', help='numeric key to order timepoints')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--cluster-colors', default=None)
    ap.add_argument('--use-raw', action='store_true', help='score from adata.raw')
    ap.add_argument('--layer', default=None, help='score from this layer instead of .X')
    ap.add_argument('--vlim', type=float, default=2.0)
    ap.add_argument('--cmap', default='RdBu_r')
    ap.add_argument('--panels', choices=['both', 'clusters', 'timepoints'], default='both')
    ap.add_argument('--orient', choices=['landscape', 'portrait'], default='landscape',
                    help='landscape = GEPs on rows, panels side-by-side (default); '
                         'portrait = GEPs on columns, groups on rows, cluster panel '
                         'stacked above timepoint panel (A4 / single-column friendly)')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig_gep_correspondence')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    import scanpy as sc
    import anndata as ad

    geps = load_geps(args.gep, args.top)
    adata = sc.read_h5ad(args.h5ad)
    print(f'[load] {adata.n_obs:,} cells x {adata.n_vars:,} genes')

    # pick the scoring matrix
    if args.layer:
        if args.layer not in adata.layers:
            raise SystemExit(f'layer {args.layer!r} not found; have {list(adata.layers)}')
        adata.X = adata.layers[args.layer]
        use_raw = False
        print(f'[score] using layer {args.layer!r}')
    elif args.use_raw:
        if adata.raw is None:
            raise SystemExit('--use-raw but adata.raw is None')
        use_raw = True
        print('[score] using adata.raw')
    else:
        use_raw = adata.raw is not None
        src = '.raw' if use_raw else '.X'
        xmin = float(adata.X.min()) if not use_raw else float(adata.raw.X.min())
        print(f'[score] using {src} (auto); min value = {xmin:.2f}')
        if not use_raw and xmin < -0.5:
            print('[warn] .X looks SCALED (negative values). Scoring still ranks groups '
                  'correctly (per-GEP z-score absorbs the scale), but for a textbook '
                  'score_genes pass re-export log-norm into .raw or pass --layer.')

    var_names = set((adata.raw.var_names if use_raw else adata.var_names))
    vn_upper = {v.upper(): v for v in var_names}
    score_cols = []
    for g, genes in geps.items():
        present = [vn_upper[x.upper()] for x in genes if x.upper() in vn_upper]
        print(f'  {g}: {len(present)}/{len(genes)} genes matched')
        if not present:
            print(f'  [skip] {g}: no genes matched')
            continue
        sc.tl.score_genes(adata, present, score_name=g, use_raw=use_raw)
        score_cols.append(g)

    obs = adata.obs
    ck, tk, wk = args.cluster_key, args.timepoint_key, args.week_key

    # ---- per-cluster matrix (GEP x cluster), z-scored per GEP across clusters ----
    cl = obs[ck].astype(str)
    per_c = pd.DataFrame({g: obs[g].groupby(cl).mean() for g in score_cols}).T  # GEP x cluster
    cats = list(per_c.columns)
    if Path(args.cluster_order).exists():
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        cols_c = [c for c in co if c in cats] + [c for c in cats if c not in co]
    else:
        cols_c = sorted(cats, key=lambda s: int(s) if s.isdigit() else 1e9)
    per_c = per_c[cols_c]
    Zc = per_c.sub(per_c.mean(axis=1), axis=0).div(per_c.std(axis=1, ddof=0).replace(0, 1), axis=0)

    # order GEP rows by peak-cluster position along the axis (diagonal layout)
    posc = {c: i for i, c in enumerate(cols_c)}
    peak_c = Zc.idxmax(axis=1)
    row_order = sorted(score_cols, key=lambda g: posc.get(peak_c[g], 99))
    Zc = Zc.loc[row_order]

    # ---- per-timepoint matrix (GEP x timepoint), z-scored per GEP across timepoints ----
    tp = obs[tk].astype(str)
    per_t = pd.DataFrame({g: obs[g].groupby(tp).mean() for g in score_cols}).T
    if wk in obs.columns:
        wmed = obs.groupby(tp)[wk].median()
        cols_t = list(wmed.sort_values().index)
    else:
        cols_t = sorted(per_t.columns)
    per_t = per_t[cols_t]
    Zt = per_t.sub(per_t.mean(axis=1), axis=0).div(per_t.std(axis=1, ddof=0).replace(0, 1), axis=0)
    Zt = Zt.loc[row_order]

    # write tables
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    per_c.loc[row_order].to_csv(outdir / f'{args.stem}_cluster_scores.tsv', sep='\t')
    per_t.loc[row_order].to_csv(outdir / f'{args.stem}_timepoint_scores.tsv', sep='\t')

    # ---- plot ----
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle

    colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if (args.cluster_colors and Path(args.cluster_colors).exists()) else {}
    norm = TwoSlopeNorm(vmin=-args.vlim, vcenter=0, vmax=args.vlim)
    ylabs = [f'{g}  {GEP_LABELS.get(g, "")}'.rstrip() for g in row_order]

    def draw(ax, Z, cols, title, mark_xcolor=False, xrot=0):
        im = ax.imshow(Z.values, aspect='auto', cmap=args.cmap, norm=norm)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, fontsize=8, rotation=xrot,
                           ha=('right' if xrot else 'center'),
                           rotation_mode='anchor')
        ax.set_yticks(range(len(row_order))); ax.set_yticklabels(ylabs, fontsize=8)
        if mark_xcolor and colors:
            for t, c in zip(ax.get_xticklabels(), cols):
                t.set_color(colors.get(c, 'black'))
        # box the per-GEP argmax group (the program's top correspondence)
        for i, g in enumerate(row_order):
            j = int(np.argmax(Z.loc[g].values))
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                   edgecolor='k', lw=1.4))
        ax.tick_params(length=0)
        ax.set_title(title, fontsize=10)
        return im

    def draw_T(ax, Z, groups, ylabel, mark_ycolor=False, show_xlabels=True):
        # transposed: GEPs on the x-axis (7 cols), groups down the rows
        im = ax.imshow(Z.values.T, aspect='auto', cmap=args.cmap, norm=norm)
        ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups, fontsize=8)
        if mark_ycolor and colors:
            for t, c in zip(ax.get_yticklabels(), groups):
                t.set_color(colors.get(c, 'black'))
        if show_xlabels:
            ax.set_xticks(range(len(row_order)))
            ax.set_xticklabels(row_order, fontsize=8, rotation=45, ha='left',
                               rotation_mode='anchor')
            ax.xaxis.set_ticks_position('top'); ax.xaxis.set_label_position('top')
        else:
            ax.set_xticks([])
        # box the per-GEP argmax group: for GEP column j, its top group (row i)
        for j, g in enumerate(row_order):
            i = int(np.argmax(Z.loc[g].values))
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                   edgecolor='k', lw=1.4))
        ax.tick_params(length=0); ax.set_ylabel(ylabel, fontsize=9)
        return im

    if args.orient == 'portrait':
        if args.panels == 'both':
            nC, nT = Zc.shape[1], Zt.shape[1]
            fig, axes = plt.subplots(
                2, 1, figsize=(0.42 * len(row_order) + 3.2, 0.30 * (nC + nT) + 2.4),
                gridspec_kw={'height_ratios': [nC, nT], 'hspace': 0.06})
            im = draw_T(axes[0], Zc, cols_c, 'cluster (adeno → NE)',
                        mark_ycolor=True, show_xlabels=True)
            draw_T(axes[1], Zt, cols_t, 'timepoint', show_xlabels=False)
            cb = fig.colorbar(im, ax=axes, shrink=0.5, aspect=14, pad=0.02)
        else:
            Z, groups, ylab, mc = (Zc, cols_c, 'cluster (adeno → NE)', True) \
                if args.panels == 'clusters' else (Zt, cols_t, 'timepoint', False)
            fig, ax = plt.subplots(figsize=(0.42 * len(row_order) + 3.2,
                                            0.30 * Z.shape[1] + 2.4))
            im = draw_T(ax, Z, groups, ylab, mark_ycolor=mc)
            cb = fig.colorbar(im, ax=ax, shrink=0.5, aspect=14, pad=0.02)
    elif args.panels == 'both':
        nC, nT = Zc.shape[1], Zt.shape[1]
        fig, axes = plt.subplots(
            1, 2, figsize=(0.34 * (nC + nT) + 4.5, 0.42 * len(row_order) + 1.8),
            gridspec_kw={'width_ratios': [nC, nT], 'wspace': 0.05})
        draw(axes[0], Zc, cols_c, 'GEP × cluster', mark_xcolor=True)
        axes[0].set_xlabel('cluster (adenocarcinoma → neuroendocrine)', fontsize=9)
        im = draw(axes[1], Zt, cols_t, 'GEP × timepoint', xrot=45)
        axes[1].set_xlabel('timepoint', fontsize=9)
        axes[1].set_yticklabels([])
        cb = fig.colorbar(im, ax=axes, shrink=0.5, aspect=14, pad=0.02)
    else:
        Z, cols, title, mx, xl = (Zc, cols_c, 'GEP × cluster', True,
                                  'cluster (adenocarcinoma → neuroendocrine)') \
            if args.panels == 'clusters' else \
            (Zt, cols_t, 'GEP × timepoint', False, 'timepoint')
        fig, ax = plt.subplots(figsize=(0.34 * Z.shape[1] + 4.0, 0.42 * len(row_order) + 1.8))
        im = draw(ax, Z, cols, title, mark_xcolor=mx, xrot=(45 if args.panels == 'timepoints' else 0))
        ax.set_xlabel(xl, fontsize=9)
        cb = fig.colorbar(im, ax=ax, shrink=0.5, aspect=14, pad=0.02)
    cb.set_label('GEP score (z across groups)', fontsize=8); cb.ax.tick_params(labelsize=7)

    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}} and *_scores.tsv to {outdir}/')


if __name__ == '__main__':
    main()
