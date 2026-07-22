#!/usr/bin/env python3
"""
plot_marker_dotplot.py — grouped marker dot plot (Fig 2C and its supplements) from a
<group>\\t<gene> TSV + the h5ad object. Self-contained renderer (no scanpy dotplot) so the
colour scaling is fully controlled.

Dot SIZE  = fraction of cells in the cluster expressing the gene (>0).
Dot COLOR = mean expression, scaled across clusters per gene:
  --scale zscore  (DEFAULT) z-score across clusters, diverging RdBu_r centred at 0,
                  limits +/- --vlim (default 2). Seurat-DotPlot convention: red = enriched
                  in that cluster, blue = depleted. Matches plot_cluster_summary.py.
  --scale var     min-max 0..1 per gene, sequential Reds.
  --scale none    raw mean lognorm, sequential Reds.

Optional --deg-table <FindAllMarkers.tsv> overlays a black '*' on every dot that is a
BH-significant POSITIVE DEG for that cluster (p_val_adj < --deg-padj and avg_log2FC > --deg-lfc).
This draws the eye to genuine cluster enrichment rather than to z-scored colour (which always
paints some cluster red — e.g. CD44 looks dark at cl17 but is NS there). Default off → Fig 2C
behaviour is unchanged.

Works for fig2c_main.tsv, fig2c_supp_canonical.tsv (groups = programs) and
fig3_datadriven.tsv (groups = cluster IDs). Genes on x grouped with the group label at the
BOTTOM (caption-matching); --swap-axes puts genes on y with group labels on the left.

CRITICAL for the LTL331 base object: .X is SCALED — fractions/means MUST come from .raw
(lognorm). Uses --use-raw by default. Run check_fig2_markers.py first.

Examples:
  python plot_marker_dotplot.py --h5ad ltl331_base.h5ad --markers data/fig2c_main.tsv \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --outdir figures --stem fig2c
  python plot_marker_dotplot.py --h5ad ltl331_base.h5ad --markers data/fig2c_supp_canonical.tsv \\
      --cluster-order data/cluster_order.json --outdir figures --stem fig2c_suppA --figwidth 22
  python plot_marker_dotplot.py --h5ad ltl331_base.h5ad --markers data/fig3_datadriven.tsv \\
      --cluster-order data/cluster_order.json --outdir figures --stem fig3 --swap-axes
"""
import argparse, json, textwrap
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def load_markers(path):
    groups = OrderedDict()
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 2:
            print(f'[warn] skipping non-2-column line: {line!r}'); continue
        g, gene = parts[0].strip(), parts[1].strip()
        groups.setdefault(g, [])
        if gene not in groups[g]:
            groups[g].append(gene)
    return groups


def compute(adata, genes, cluster_key, order, use_raw):
    """mean (lognorm) and fraction-expressing per (cluster, gene). rows=order, cols=present."""
    src = adata.raw if (use_raw and adata.raw is not None) else adata
    var = list(map(str, src.var_names)); idx = {g: i for i, g in enumerate(var)}
    present = [g for g in genes if g in idx]
    missing = [g for g in genes if g not in idx]
    if not present:
        return pd.DataFrame(index=order), pd.DataFrame(index=order), missing
    gi = [idx[g] for g in present]
    X = src.X[:, gi]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    cl = adata.obs[cluster_key].astype(str).values
    means, fracs = [], []
    for c in order:
        m = cl == c
        if m.sum() == 0:
            means.append([np.nan] * len(present)); fracs.append([np.nan] * len(present)); continue
        sub = X[m]
        means.append(sub.mean(axis=0)); fracs.append((sub > 0).mean(axis=0))
    return (pd.DataFrame(means, index=order, columns=present),
            pd.DataFrame(fracs, index=order, columns=present), missing)


def wrap_group_label(lab, width_cols):
    """Wrap a group descriptor so it fits a narrow gene span (prevents adjacent-label overlap
    when a group has few genes). A literal '|' is a manual line break and takes precedence;
    otherwise auto-wrap at ~3 chars per gene column."""
    if '|' in lab:
        return lab.replace('|', '\n')
    maxchars = max(5, int(width_cols * 3))
    if len(lab) <= maxchars:
        return lab
    return textwrap.fill(lab, maxchars, break_long_words=False, break_on_hyphens=False)


def load_sig_degs(path, padj_thr, lfc_thr):
    """Set of (cluster_str, gene) that are BH-significant positive DEGs in a FindAllMarkers TSV
    (cols: avg_log2FC, p_val_adj, cluster, gene). Tolerant of CRLF / stray whitespace."""
    df = pd.read_csv(path, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    for c in ('gene', 'cluster'):
        df[c] = df[c].astype(str).str.strip()
    sig = df[(df['p_val_adj'] < padj_thr) & (df['avg_log2FC'] > lfc_thr)]
    return set(zip(sig['cluster'], sig['gene']))


def scale_cols(mean_df, how):
    if how == 'zscore':
        mu = mean_df.mean(axis=0); sd = mean_df.std(axis=0, ddof=0).replace(0, 1)
        return (mean_df - mu) / sd
    if how == 'var':
        lo = mean_df.min(axis=0); rng = (mean_df.max(axis=0) - lo).replace(0, 1)
        return (mean_df - lo) / rng
    return mean_df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--markers', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--cluster-order', default=None)
    ap.add_argument('--reverse-clusters', action='store_true',
                    help='reverse the cluster order (e.g. NE->adeno instead of adeno->NE). '
                         'Flips the whole axis consistently; works with or without --cluster-order.')
    ap.add_argument('--cluster-colors', default=None,
                    help='cluster->hex json (e.g. data/cluster_colors_18.json) to colour the '
                         'cluster tick labels')
    ap.add_argument('--use-raw', dest='use_raw', action='store_true', default=True)
    ap.add_argument('--no-use-raw', dest='use_raw', action='store_false')
    ap.add_argument('--scale', choices=['zscore', 'var', 'none'], default='zscore')
    ap.add_argument('--vlim', type=float, default=2.0, help='z-score colour limit (default 2)')
    ap.add_argument('--deg-table', default=None,
                    help='FindAllMarkers TSV; overlay a * on BH-significant positive DEGs '
                         '(p_val_adj<--deg-padj, avg_log2FC>--deg-lfc). Default off.')
    ap.add_argument('--deg-padj', type=float, default=0.05)
    ap.add_argument('--deg-lfc', type=float, default=0.0)
    ap.add_argument('--group-label-rotation', type=float, default=0.0,
                    help='rotate group descriptors (deg); 0 = horizontal (default). Labels are '
                         'auto-wrapped to fit narrow gene spans; use "|" in a TSV group name to '
                         'force a line break.')
    ap.add_argument('--group-label-size', type=float, default=8.0)
    ap.add_argument('--swap-axes', action='store_true', help='genes on y, clusters on x')
    ap.add_argument('--smin', type=float, default=4.0)
    ap.add_argument('--smax', type=float, default=170.0)
    ap.add_argument('--figwidth', type=float, default=None)
    ap.add_argument('--figheight', type=float, default=None)
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='marker_dotplot')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm, Normalize
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    cats = sorted(adata.obs[args.cluster_key].unique(),
                  key=lambda s: int(s) if s.isdigit() else 1e9)
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        order = [c for c in co if c in cats] + [c for c in cats if c not in co]
    else:
        order = cats
    if args.reverse_clusters:
        order = order[::-1]

    sig_set = set()
    if args.deg_table:
        sig_set = load_sig_degs(args.deg_table, args.deg_padj, args.deg_lfc)
        print(f'[info] {len(sig_set)} significant (cluster,gene) DEGs loaded for * overlay '
              f'(padj<{args.deg_padj}, lfc>{args.deg_lfc})')

    groups = load_markers(args.markers)
    # flatten to gene list + group spans, dropping genes absent from the object
    flat, gene_group, missing_all = [], [], []
    for grp, genes in groups.items():
        m, f, missing = compute(adata, genes, args.cluster_key, order, args.use_raw)
        missing_all += missing
        for g in m.columns:
            flat.append(g); gene_group.append(grp)
    if missing_all:
        print(f'[warn] {len(missing_all)} symbol(s) absent from object and DROPPED: {missing_all}')
        print('       -> run check_fig2_markers.py for the correct symbols.')
    if not flat:
        raise SystemExit('no marker genes resolve in the object')

    mean_all, frac_all, _ = compute(adata, flat, args.cluster_key, order, args.use_raw)
    mean_all = mean_all[flat]; frac_all = frac_all[flat]
    color_df = scale_cols(mean_all, args.scale)

    # colour scheme
    if args.scale == 'zscore':
        cmap = plt.get_cmap('RdBu_r'); norm = TwoSlopeNorm(vmin=-args.vlim, vcenter=0, vmax=args.vlim)
        cbar_label = 'scaled mean expression (z-score)'
    elif args.scale == 'var':
        cmap = plt.get_cmap('Reds'); norm = Normalize(0, 1)
        cbar_label = 'relative expression (0–1 per gene)'
    else:
        cmap = plt.get_cmap('Reds'); norm = Normalize(float(np.nanmin(color_df.values)),
                                                      float(np.nanmax(color_df.values)))
        cbar_label = 'mean expression (lognorm)'

    cl_colors = {}
    if args.cluster_colors:
        cl_colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()}

    ng, nc = len(flat), len(order)
    genes_on_x = not args.swap_axes

    # ---- figure ----
    if genes_on_x:
        W = args.figwidth or (0.23 * ng + 3.5)
        H = args.figheight or (0.30 * nc + 2.5)
    else:
        W = args.figwidth or (0.30 * nc + 3.5)
        H = args.figheight or (0.23 * ng + 2.5)
    fig, ax = plt.subplots(figsize=(W, H))

    def sz(frac):
        return args.smin + np.clip(frac, 0, 1) * (args.smax - args.smin)

    # group boundaries (cumulative gene counts per group, in flat order)
    bounds, labels_g, c = [], [], 0
    cur = gene_group[0]; start = 0
    for i, gg in enumerate(gene_group + [None]):
        if gg != cur:
            bounds.append((start, i - 1)); labels_g.append(cur); start = i; cur = gg

    for gi_, g in enumerate(flat):
        for ci_, cl in enumerate(order):
            v = color_df.iloc[ci_][g]; f = frac_all.iloc[ci_][g]
            if pd.isna(v):
                continue
            if genes_on_x:
                x, y = gi_, nc - 1 - ci_
            else:
                x, y = ci_, ng - 1 - gi_
            ax.scatter(x, y, s=sz(f), color=cmap(norm(v)), edgecolor='black',
                       lw=0.3, zorder=3)
            if sig_set and (cl, g) in sig_set:
                ax.scatter(x, y, marker='*', s=34, c='black', edgecolors='white',
                           linewidths=0.4, zorder=5)

    if genes_on_x:
        ax.set_xlim(-0.5, ng - 0.5); ax.set_ylim(-0.5, nc - 0.5)
        ax.set_xticks(range(ng)); ax.set_xticklabels(flat, rotation=90, fontsize=12)#6.5)
        ax.set_yticks(range(nc)); ax.set_yticklabels(list(reversed(order)), fontsize=12)#7)
        if cl_colors:
            for t, cl in zip(ax.get_yticklabels(), reversed(order)):
                t.set_color(cl_colors.get(cl, 'black'))
        # group separators + labels along the bottom
        tr = ax.get_xaxis_transform()
        for (s, e) in bounds[:-1]:
            ax.axvline(e + 0.5, color='0.6', lw=0.8, zorder=1)
        for (s, e), lab in zip(bounds, labels_g):
            ax.text((s + e) / 2, -0.16, wrap_group_label(lab, e - s + 1), transform=tr,
                    ha='center', va='top', fontsize=args.group_label_size, fontweight='bold',
                    rotation=args.group_label_rotation,
                    rotation_mode='anchor' if args.group_label_rotation else 'default')
    else:
        ax.set_xlim(-0.5, nc - 0.5); ax.set_ylim(-0.5, ng - 0.5)
        ax.set_xticks(range(nc)); ax.set_xticklabels(order, rotation=90, fontsize=7)
        if cl_colors:
            for t, cl in zip(ax.get_xticklabels(), order):
                t.set_color(cl_colors.get(cl, 'black'))
        ax.set_yticks(range(ng)); ax.set_yticklabels(list(reversed(flat)), fontsize=7.5)
        tr = ax.get_yaxis_transform()
        for (s, e) in bounds[:-1]:
            ax.axhline(ng - 1 - (e + 0.5), color='0.6', lw=0.8, zorder=1)
        for (s, e), lab in zip(bounds, labels_g):
            ax.text(-0.16, ng - 1 - (s + e) / 2, wrap_group_label(lab, e - s + 1), transform=tr,
                    ha='right', va='center', fontsize=args.group_label_size, fontweight='bold',
                    rotation=90)

    ax.tick_params(length=0)
    for sp_ in ('top', 'right', 'left', 'bottom'):
        ax.spines[sp_].set_visible(False)

    # colorbar — fraction=0.035 (was default 0.15, which stole 15% of the axis width) gives
    # that width back to the dot plot; label kept LARGE for readability.
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    # thinner bar (aspect 16) narrows the legend column; label kept LARGE (13 pt) for readability.
    cb = fig.colorbar(sm, ax=ax, fraction=0.022, shrink=0.5, aspect=16, pad=0.01)
    cb.set_label(cbar_label, fontsize=13); cb.ax.tick_params(labelsize=8)

    # size legend — narrow footprint (tight spacing/pad)
    for frac in (0.25, 0.5, 0.75, 1.0):
        ax.scatter([], [], s=sz(frac), color='0.4', edgecolor='black', lw=0.3,
                   label=f'{int(frac*100)}%')
    leg_size = ax.legend(title='% expressing', loc='upper left', bbox_to_anchor=(1.005, 1.0),
                         fontsize=7.5, title_fontsize=9, frameon=False, labelspacing=0.5,
                         borderpad=0.25, handletextpad=0.3)
    if sig_set:
        star = ax.scatter([], [], marker='*', s=34, c='black', edgecolors='white',
                          linewidths=0.4, label=f'BH-sig DEG (padj<{args.deg_padj:g})')
        ax.legend(handles=[star], loc='lower left', bbox_to_anchor=(1.005, 0.0),
                  fontsize=7.5, frameon=False, borderpad=0.3)
        ax.add_artist(leg_size)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}} to {outdir}/'
          f'  ({ng} genes x {nc} clusters, scale={args.scale})')


if __name__ == '__main__':
    main()
