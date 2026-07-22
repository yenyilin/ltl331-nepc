#!/usr/bin/env python3
"""
plot_origin_vs_harmony.py

Track where original clusters land in the joint Harmony clustering
for each dataset as a single stacked figure:

  [STAGE RIBBON] (optional, --cluster-order)
               joint clusters colored/labelled by the phenotypic stage
               (AR-high PRAD -> ... -> NE-terminal) of their dominant LTL331
               original cluster. Makes the x-axis read as the adeno->NE axis.

  TOP panel    rows = dataset_originalCluster (LTL331, then Gao, then Li),
               cols = Harmony joint cluster,
               colour = row-normalised fraction (each original cluster's row
               sums to 1 -> "where does this cluster go").

  COHORT panel 3 rows = cohort (LTL331 / Gao / Li),
               cols = same Harmony clusters,
               colour + annotation = % of each Harmony cluster's cells from
               that cohort (each column sums to 100 -> cross-cohort make-up).

  [PHENOTYPE panel] (optional, --phenotype-col)
               rows = clinical disease phenotype (Gao/Li only),
               cols = same Harmony clusters,
               colour = % of each cluster's CLINICAL cells from that phenotype.
               Shows which clinical disease states populate each joint cluster.

All panels share the x-axis. With --cluster-order, columns and rows are ordered
along the adenocarcinoma->neuroendocrine phenotypic axis (data/cluster_order.json),
so the cross-cohort co-clustering reads as a clinical progression.

Usage:
    python plot_origin_vs_harmony.py \\
        --h5ad integrated.h5ad --outdir origin_map \\
        --res-col RNA_snn_res.0.55 --group-col dataset \\
        --original-col original_cluster --ltl331-label LTL331 \\
        --cluster-order data/cluster_order.json \\
        --phenotype-col patient_phenotype \\
        --plot-format pdf png
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from _utils import setup_style, numkey as numeric_key

# ColorBrewer Set2 (categorical, colorblind-safe)
COHORT_PALETTE = {'LTL331': '#66C2A5', 'Gao': '#FC8D62', 'Li': '#8DA0CB'}

# Stage ribbon colors along the adeno->NE axis (Intermediate = cl17 pink, matching
# cluster_colors_18); 'clinical-only' = grey for joint clusters with no LTL331 anchor.
STAGE_PALETTE = {
    'AR-high PRAD': '#2166ac', 'AR-low PRAD': '#67a9cf', 'pre-EMT': '#d1e5f0',
    'Intermediate': '#e7298a', 'NE-early': '#fdae61', 'NE-mid': '#f46d43',
    'NE-terminal': '#a50026', 'clinical-only': '#cccccc',
}


def load_cluster_order(path):
    """Return (pos: cluster->axis index, stage_of: cluster->stage, axis_label)."""
    d = json.load(open(path))
    order = [str(c) for c in d['order']]
    pos = {c: i for i, c in enumerate(order)}
    stage_of = {}
    for stage, members in d.get('stages', {}).items():
        for c in members:
            stage_of[str(c)] = stage
    return pos, stage_of, d.get('axis_label', 'cluster (adenocarcinoma → neuroendocrine)')


def ltl_axis_maps(obs, group_col, orig_col, res_col, ltl331_label, pos):
    """Per joint cluster: mean LTL331 axis position (for ordering) and dominant
    LTL331 original cluster (for the stage ribbon)."""
    ltl = obs[obs[group_col].astype(str) == ltl331_label].copy()
    if ltl.empty:
        return {}, {}
    ltl[res_col] = ltl[res_col].astype(str)
    ltl[orig_col] = ltl[orig_col].astype(str)
    ltl['_pos'] = ltl[orig_col].map(pos)
    mean_pos = ltl.groupby(res_col)['_pos'].mean().to_dict()       # joint -> mean axis pos
    ct = pd.crosstab(ltl[res_col], ltl[orig_col])
    dom = ct.idxmax(axis=1).to_dict()                              # joint -> dominant LTL orig
    return mean_pos, dom


def build_matrices(obs, group_col, orig_col, res_col, ltl331_label, pheno_col=None):
    """Return (row-frac origin x harmony, cohort % comp, clinical-phenotype % or None)."""
    obs = obs.copy()
    for c in (group_col, orig_col, res_col):
        obs[c] = obs[c].astype(str)
    obs['orig_label'] = obs[group_col] + '_' + obs[orig_col]

    ct = pd.crosstab(obs['orig_label'], obs[res_col])
    frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)

    coh = pd.crosstab(obs[group_col], obs[res_col])
    comp = coh.div(coh.sum(axis=0).replace(0, np.nan), axis=1) * 100

    pheno = None
    if pheno_col is not None:
        clin = obs[obs[group_col] != ltl331_label].copy()         # clinical cells only
        clin[pheno_col] = clin[pheno_col].astype(str)
        pc = pd.crosstab(clin[pheno_col], clin[res_col])
        pheno = pc.div(pc.sum(axis=0).replace(0, np.nan), axis=1) * 100

    # default ordering: harmony cols numeric; origin rows by cohort then numeric
    cols = sorted(frac.columns, key=numeric_key)
    frac = frac[cols]
    comp = comp.reindex(columns=cols)
    if pheno is not None:
        pheno = pheno.reindex(columns=cols)
    ds_rank = {ltl331_label: 0, 'LTL331': 0, 'other': 0, 'Gao': 1, 'Li': 2}

    def row_key(lbl):
        ds, cl = lbl.rsplit('_', 1)
        return (ds_rank.get(ds, 99), numeric_key(cl))
    frac = frac.loc[sorted(frac.index, key=row_key)]
    comp = comp.loc[sorted(comp.index, key=lambda d: ds_rank.get(d, 99))]
    comp = comp.rename(index={ltl331_label: 'LTL331', 'other': 'LTL331'})
    return frac, comp, pheno


def apply_phenotypic_order(frac, comp, pheno, obs, group_col, orig_col, res_col,
                           ltl331_label, pos, stage_of):
    """Reorder columns + rows along the adeno->NE axis; return ordering + per-column
    stage list for the ribbon."""
    mean_pos, dom = ltl_axis_maps(obs, group_col, orig_col, res_col, ltl331_label, pos)
    big = (max(pos.values()) + 2) if pos else 0

    # columns: by mean LTL331 axis pos; clinical-only clusters (no LTL cells) -> NE end
    col_order = sorted(frac.columns,
                       key=lambda c: (mean_pos.get(c, big + 1), numeric_key(c)))
    frac = frac[col_order]
    comp = comp.reindex(columns=col_order)
    if pheno is not None:
        pheno = pheno.reindex(columns=col_order)
    col_ord = {c: i for i, c in enumerate(col_order)}

    # rows of frac: keep cohort blocks; within block, by axis position
    ds_rank = {ltl331_label: 0, 'LTL331': 0, 'other': 0, 'Gao': 1, 'Li': 2}

    def landing(lbl):
        ds, cl = lbl.rsplit('_', 1)
        if ds_rank.get(ds, 99) == 0 and cl in pos:        # LTL331 row -> its own axis pos
            return (0, pos[cl])
        row = frac.loc[lbl].fillna(0).values              # clinical -> weighted landing pos
        w = row.sum()
        lp = float(np.dot(row, np.arange(len(row))) / w) if w > 0 else big
        return (ds_rank.get(ds, 99), lp)
    frac = frac.loc[sorted(frac.index, key=landing)]

    # phenotype rows: by weighted-mean column position (data-driven disease axis)
    if pheno is not None and pheno.shape[0] > 1:
        def pheno_key(r):
            row = pheno.loc[r].fillna(0).values
            w = row.sum()
            return float(np.dot(row, np.arange(len(row))) / w) if w > 0 else 1e9
        pheno = pheno.loc[sorted(pheno.index, key=pheno_key)]

    stage_per_col = [stage_of.get(str(dom.get(c)), 'clinical-only') for c in col_order]
    return frac, comp, pheno, stage_per_col


def _resolve_cmaps():
    import matplotlib.pyplot as plt
    return plt.get_cmap('YlGnBu'), plt.get_cmap('YlOrBr')


def plot_stacked(frac, comp, outdir, prefix, formats, pheno=None,
                 stage_per_col=None, axis_label='Harmony joint cluster'):
    plt = setup_style(font_size=7)
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Patch

    n_orig, n_cols = frac.shape
    n_coh = comp.shape[0]
    n_ph = pheno.shape[0] if pheno is not None else 0
    has_ribbon = stage_per_col is not None
    row_ds = [lbl.rsplit('_', 1)[0] for lbl in frac.index]

    nz = frac.values[frac.values > 0]
    top_vmax = max(float(np.quantile(nz, 0.95)) if nz.size else 1.0, 0.05)

    fig_w = max(6.0, n_cols * 0.40 + 4.0)
    fig_h = max(6.5, n_orig * 0.16 + n_coh * 0.55 + n_ph * 0.30 + 2.8 + 0.3 * has_ribbon)

    # dynamic gridspec: [cbar_top][ribbon?][top][cbar_bot][comp][cbar_ph?][pheno?]
    heights, slots = [], []
    heights.append(0.22); slots.append('cbar_top')
    if has_ribbon:
        heights.append(0.30); slots.append('ribbon')
    heights.append(n_orig); slots.append('top')
    heights.append(0.22); slots.append('cbar_bot')
    heights.append(max(n_coh * 1.4, 4)); slots.append('comp')
    if pheno is not None:
        heights.append(0.22); slots.append('cbar_ph')
        heights.append(max(n_ph * 1.2, 3)); slots.append('pheno')

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=0.12)
    ax = {}
    for i, name in enumerate(slots):
        ax[name] = fig.add_subplot(gs[i, 0])
    # share x with the top heatmap
    anchor = ax['top']
    for name in ('ribbon', 'comp', 'pheno'):
        if name in ax:
            ax[name].sharex(anchor)

    cmap_top, cmap_bot = _resolve_cmaps()
    cmap_top = cmap_top.copy(); cmap_top.set_bad('#ededed')

    # ----- optional stage ribbon ------------------------------------------- #
    if has_ribbon:
        idx = {s: k for k, s in enumerate(STAGE_PALETTE)}
        from matplotlib.colors import ListedColormap, BoundaryNorm
        ribbon_cmap = ListedColormap([STAGE_PALETTE[s] for s in STAGE_PALETTE])
        data = np.array([[idx[s] for s in stage_per_col]])
        ax['ribbon'].imshow(data, aspect='auto', cmap=ribbon_cmap,
                            norm=BoundaryNorm(range(len(STAGE_PALETTE) + 1), len(STAGE_PALETTE)),
                            interpolation='nearest', rasterized=True)
        # label contiguous stage runs centered over their span
        s0 = 0
        for j in range(1, n_cols + 1):
            if j == n_cols or stage_per_col[j] != stage_per_col[s0]:
                lab = stage_per_col[s0]
                ax['ribbon'].text((s0 + j - 1) / 2, 0, lab, ha='center', va='center',
                                  fontsize=5, rotation=0,
                                  color='white' if lab in ('AR-high PRAD', 'NE-terminal',
                                                            'Intermediate') else '#222')
                s0 = j
        ax['ribbon'].set_yticks([]); ax['ribbon'].tick_params(labelbottom=False, length=0)
        ax['ribbon'].set_ylabel('stage', fontsize=6, rotation=0, ha='right', va='center')
        for s in ax['ribbon'].spines.values():
            s.set_visible(False)

    # ----- top panel: origin x harmony row fraction ------------------------ #
    top_data = np.ma.masked_equal(frac.values, 0.0)
    im1 = ax['top'].imshow(top_data, aspect='auto', cmap=cmap_top,
                           norm=PowerNorm(gamma=0.5, vmin=1e-6, vmax=top_vmax),
                           interpolation='nearest')
    ax['top'].set_yticks(range(n_orig))
    ax['top'].set_yticklabels(frac.index, fontsize=5)
    for tick, ds in zip(ax['top'].get_yticklabels(), row_ds):
        tick.set_color(COHORT_PALETTE.get(ds, '#444')); tick.set_fontweight('bold')
    ax['top'].set_ylabel('dataset_originalCluster')
    ax['top'].tick_params(labelbottom=False, length=0)
    for s in ax['top'].spines.values():
        s.set_visible(False)
    prev = None
    for i, ds in enumerate(row_ds):
        if prev is not None and ds != prev:
            ax['top'].axhline(i - 0.5, color='#333333', lw=0.6)
        prev = ds
    cb1 = fig.colorbar(im1, cax=ax['cbar_top'], orientation='horizontal')
    cb1.set_label(f'frac of original cluster (√-scaled, clipped at 95th pct = '
                  f'{top_vmax:.2f}; 0 → grey)', fontsize=6)
    cb1.ax.tick_params(labelsize=6, length=0)
    ax['cbar_top'].xaxis.set_ticks_position('top'); ax['cbar_top'].xaxis.set_label_position('top')

    # ----- cohort composition panel ---------------------------------------- #
    im2 = ax['comp'].imshow(comp.values, aspect='auto', cmap=cmap_bot,
                            vmin=0, vmax=100, interpolation='nearest')
    ax['comp'].set_yticks(range(n_coh)); ax['comp'].set_yticklabels(comp.index)
    for tick, coh in zip(ax['comp'].get_yticklabels(), comp.index):
        tick.set_color(COHORT_PALETTE.get(coh, '#444')); tick.set_fontweight('bold')
    bottom_ax = ax['pheno'] if pheno is not None else ax['comp']
    ax['comp'].set_ylabel('cohort %')
    ax['comp'].tick_params(length=0)
    if pheno is None:
        xtick_rot = 0 if n_cols <= 18 else 90
        ax['comp'].set_xticks(range(n_cols))
        ax['comp'].set_xticklabels(comp.columns, rotation=xtick_rot, fontsize=6)
        ax['comp'].set_xlabel(axis_label)
    else:
        ax['comp'].tick_params(labelbottom=False)
    for s in ax['comp'].spines.values():
        s.set_visible(False)
    cb2 = fig.colorbar(im2, cax=ax['cbar_bot'], orientation='horizontal')
    cb2.set_label('% of cluster', fontsize=6); cb2.ax.tick_params(labelsize=6, length=0)
    ax['cbar_bot'].xaxis.set_ticks_position('top'); ax['cbar_bot'].xaxis.set_label_position('top')
    for i in range(comp.shape[0]):
        for j in range(comp.shape[1]):
            v = comp.values[i, j]
            if np.isfinite(v):
                ax['comp'].text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=5.5,
                                color='white' if v >= 60 else '#3a3a3a')

    # ----- optional clinical-phenotype panel ------------------------------- #
    if pheno is not None:
        im3 = ax['pheno'].imshow(pheno.values, aspect='auto', cmap='BuPu',
                                 vmin=0, vmax=100, interpolation='nearest')
        ax['pheno'].set_yticks(range(n_ph)); ax['pheno'].set_yticklabels(pheno.index, fontsize=5.5)
        xtick_rot = 0 if n_cols <= 18 else 90
        ax['pheno'].set_xticks(range(n_cols))
        ax['pheno'].set_xticklabels(pheno.columns, rotation=xtick_rot, fontsize=6)
        ax['pheno'].set_xlabel(axis_label)
        ax['pheno'].set_ylabel('clinical\nphenotype %', fontsize=6)
        ax['pheno'].tick_params(length=0)
        for s in ax['pheno'].spines.values():
            s.set_visible(False)
        cb3 = fig.colorbar(im3, cax=ax['cbar_ph'], orientation='horizontal')
        cb3.set_label('% of cluster\'s clinical cells', fontsize=6)
        cb3.ax.tick_params(labelsize=6, length=0)
        ax['cbar_ph'].xaxis.set_ticks_position('top'); ax['cbar_ph'].xaxis.set_label_position('top')
        for i in range(pheno.shape[0]):
            for j in range(pheno.shape[1]):
                v = pheno.values[i, j]
                if np.isfinite(v) and v >= 1:
                    ax['pheno'].text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=5,
                                     color='white' if v >= 60 else '#3a3a3a')

    for fmt in formats:
        path = Path(outdir) / f'{prefix}.{fmt}'
        fig.savefig(path, format=fmt, bbox_inches='tight')
        print(f'  wrote {path}')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=__doc__)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--outdir', default='origin_map')
    p.add_argument('--res-col', default='RNA_snn_res.0.55',
                   help='obs column with the joint Harmony clustering')
    p.add_argument('--group-col', default='dataset')
    p.add_argument('--original-col', default='original_cluster')
    p.add_argument('--ltl331-label', default='LTL331',
                   help='value in --group-col identifying LTL331 cells')
    p.add_argument('--cluster-order', default=None,
                   help='path to cluster_order.json; orders columns/rows along the '
                        'adeno->NE phenotypic axis and adds a stage ribbon')
    p.add_argument('--phenotype-col', default=None,
                   help='obs column with clinical patient phenotype (e.g. '
                        'patient_phenotype, ned_status); adds a clinical-phenotype panel')
    p.add_argument('--prefix', default='origin_vs_harmony')
    p.add_argument('--plot-format', nargs='+', default=['pdf'], choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    need = [args.res_col, args.group_col, args.original_col]
    if args.phenotype_col:
        need.append(args.phenotype_col)
    for col in need:
        if col not in adata.obs.columns:
            raise ValueError(f'column "{col}" not in obs (have: {list(adata.obs.columns)})')

    frac, comp, pheno = build_matrices(
        adata.obs, args.group_col, args.original_col, args.res_col,
        args.ltl331_label, pheno_col=args.phenotype_col)

    axis_label = 'Harmony joint cluster'
    stage_per_col = None
    if args.cluster_order:
        pos, stage_of, axis_label = load_cluster_order(args.cluster_order)
        frac, comp, pheno, stage_per_col = apply_phenotypic_order(
            frac, comp, pheno, adata.obs, args.group_col, args.original_col,
            args.res_col, args.ltl331_label, pos, stage_of)
        axis_label = f'Harmony joint cluster  ({axis_label})'

    frac.to_csv(outdir / f'{args.prefix}_rowfrac.tsv', sep='\t')
    comp.to_csv(outdir / f'{args.prefix}_cohort_composition.tsv', sep='\t')
    if pheno is not None:
        pheno.to_csv(outdir / f'{args.prefix}_clinical_phenotype.tsv', sep='\t')
    print(f'Saved matrices to {outdir}/ ({frac.shape[0]} original x {frac.shape[1]} '
          f'harmony; {comp.shape[0]} cohorts'
          + (f'; {pheno.shape[0]} phenotypes' if pheno is not None else '') + ')')

    print(f'Rendering figure ({", ".join(args.plot_format)}):')
    plot_stacked(frac, comp, outdir, args.prefix, args.plot_format,
                 pheno=pheno, stage_per_col=stage_per_col, axis_label=axis_label)
    print('Done.')


if __name__ == '__main__':
    main()
