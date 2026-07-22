#!/usr/bin/env python3
"""
plot_pseudotime_qc.py — pseudotime QC supplementary figure.

In our LTL331/R longitudinal data, a natural question about any pseudotime ordering
is how it relates to real (sampling) time. Velocity pseudotime here places some wk8-12
cells at lower values than pre-castration cells, which looks wrong if pseudotime is read
as a clock. It isn't: pseudotime orders cells by transcriptional progression along the
trajectory, not by the week a sample was collected, so the two axes can disagree by design.
This script makes the pseudotime-vs-real-time dissociation explicit at the per-cell level,
so the ordering is interpreted correctly.

Panels
  A — UMAP coloured by velocity_pseudotime (sanity check: pt continuous on
      embedding, root at PRAD, tip at NE)
  B — pt-vs-real-time per-cell scatter; x = velocity_pseudotime,
      y = numeric sample week (with y-jitter), points coloured by cluster.
      The off-diagonal cells (high real-time, low pt — e.g. cluster 0 at
      wk8-12 with pt ~0.02) become visually obvious.
  C — (optional, if --latent-time-key is present) pt-vs-latent-time per-cell
      scatter coloured by cluster; reports the global Spearman ρ in the
      title so the direction-validation result is on the figure.
  D — (optional, if --cluster-order is present) pseudotime distribution per
      cluster drawn in the fixed adenocarcinoma->neuroendocrine order, with the
      Jonckheere-Terpstra ordered-trend statistics (per-cell PI + Spearman/
      Kendall, and the independence-respecting timepoint-median test) annotated.
      This panel visualises the monotone-trend result, pseudotime increasing
      along the adenocarcinoma-to-neuroendocrine cluster order; the reported
      statistics match jonckheere_trend.py. When A-D are all present the figure
      is laid out as a 2x2 grid (A B / C D).


Outputs: figures/pseudotime_qc.{pdf,png}. PDF text fonttype 42 (editable).
"""
import argparse
import json
import textwrap
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, kendalltau, mannwhitneyu, norm


# Cluster → trajectory-phase map (must match plot_cluster_summary.py)
PHASE_MAP = {
    '11': 'preCx', '15': 'preCx',
    '2':  'wk2-4', '5':  'wk2-4', '8':  'wk2-4',
    '0':  'wk8-12', '16': 'wk8-12', '12': 'wk8-12',
    '4':  'wk16',  '6':  'wk16',
    '17': 'bridge',
    '13': 'NE_entry', '14': 'NE_entry',
    '10': 'NE_ASCL1+', '1': 'NE_ASCL1-', '3': 'NE_ASCL1-',
    '9':  'NE_ASCL1-', '7': 'NE_ASCL1-',
}
# Phase super-groups for the per-phase ρ table
PHASE_GROUPS = {
    'PRAD (preCx+regression)': {'preCx', 'wk2-4', 'wk8-12', 'wk16'},
    'bridge (17)':              {'bridge'},
    'NE block':                 {'NE_entry', 'NE_ASCL1+', 'NE_ASCL1-'},
}


def parse_week(s):
    """preCx → 0; week02 / wk2 / 2 → 2; week20R → 20 (treats regrowth as the
    same numeric week)."""
    s = str(s)
    if 'preCx' in s or 'preCx1' in s.lower():
        return 0
    digits = ''.join(c for c in s if c.isdigit())
    return int(digits) if digits else np.nan


def cluster_palette(clusters):
    """Stable cluster → colour mapping (tab20 cycled if > 20 clusters)."""
    base = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors)
    return {cl: base[i % len(base)] for i, cl in enumerate(clusters)}


def subsample(adata, cap, rng):
    if cap is None or adata.n_obs <= cap:
        return np.arange(adata.n_obs)
    return rng.choice(adata.n_obs, cap, replace=False)


# --------------------------------------------------------------------------- #
# Jonckheere-Terpstra ordered-trend test (copied verbatim from
# jonckheere_trend.py so panel D reports the same PI/Z/p as the stats script).
# --------------------------------------------------------------------------- #
def jt_variance_tiecorrected(group_sizes, tie_sizes, N):
    """Exact JT null variance with tie correction (Hollander, Wolfe & Chicken)."""
    n = np.asarray(group_sizes, dtype=float)
    e = np.asarray(tie_sizes, dtype=float)
    if N < 3:
        return np.nan
    term1 = (N * (N - 1) * (2 * N + 5)
             - np.sum(n * (n - 1) * (2 * n + 5))
             - np.sum(e * (e - 1) * (2 * e + 5))) / 72.0
    term2 = (np.sum(n * (n - 1) * (n - 2)) * np.sum(e * (e - 1) * (e - 2))) \
        / (36.0 * N * (N - 1) * (N - 2))
    term3 = (np.sum(n * (n - 1)) * np.sum(e * (e - 1))) / (8.0 * N * (N - 1))
    return term1 + term2 + term3


def jonckheere(groups):
    """groups: list of 1D arrays in the hypothesized INCREASING order.
    Returns J, tie-corrected Z, one-sided p (increasing), probabilistic index PI."""
    k = len(groups)
    J = 0.0
    for i, j in combinations(range(k), 2):
        J += mannwhitneyu(groups[j], groups[i], method="asymptotic").statistic
    n = np.array([len(g) for g in groups], dtype=float)
    N = int(n.sum())
    EJ = (N ** 2 - np.sum(n ** 2)) / 4.0
    _, tie_counts = np.unique(np.concatenate(groups), return_counts=True)
    VJ = jt_variance_tiecorrected(n, tie_counts, N)
    Z = (J - EJ) / np.sqrt(VJ) if VJ and VJ > 0 else np.nan
    denom = float(np.sum([n[a] * n[b] for a, b in combinations(range(k), 2)]))
    PI = J / denom if denom else np.nan
    P = norm.sf(Z) if np.isfinite(Z) else np.nan
    return dict(J=J, EJ=EJ, VJ=VJ, Z=Z, P=P, PI=PI, N=N, k=k)


def load_cluster_order(path):
    """Return the canonical AR->NE cluster order as a list of str cluster ids."""
    with open(path) as fh:
        obj = json.load(fh)
    if isinstance(obj, dict):
        obj = obj["order"] if "order" in obj else \
            [k for k, _ in sorted(obj.items(), key=lambda kv: kv[1])]
    return [str(x) for x in obj]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--cluster-key', default='seurat_clusters')
    p.add_argument('--pt-key', default='velocity_pseudotime')
    p.add_argument('--timepoint-key', default='sample')
    p.add_argument('--week-num-key', default='week_num',
                   help='obs column with a numeric week (preferred over '
                        'parsing the timepoint string). If absent, parse_week '
                        'is applied to --timepoint-key as a fallback.')
    p.add_argument('--latent-time-key', default='latent_time',
                   help='Optional: scVelo latent_time column for panel C; '
                        'omit panel if absent.')
    p.add_argument('--cluster-order', default='data/cluster_order.json',
                   help='JSON with the canonical AR->NE cluster order for panel D '
                        '(the Jonckheere-Terpstra monotone-trend panel). Omit the '
                        'panel if the file is absent.')
    p.add_argument('--exclude-clusters', nargs='*', default=['18'],
                   help='clusters excluded from the panel-D trend test '
                        '(default 18 = inferCNV reference, not a tumor cluster).')
    p.add_argument('--palette', default='data/cluster_colors_18.json',
                   help='JSON mapping cluster id -> hex colour (the manuscript-wide '
                        'cluster_colors_18.json). Clusters absent from the file '
                        '(or if the file is missing) fall back to the tab20 default.')
    p.add_argument('--umap-key', default='X_umap',
                   help='obsm key for UMAP coordinates')
    p.add_argument('--subsample', type=int, default=15000,
                   help='Max cells plotted per scatter panel (default 15000)')
    p.add_argument('--outdir', default='figures')
    p.add_argument('--stem', default='pseudotime_qc')
    p.add_argument('--figwidth-mm', type=float, default=None,
                   help='author the figure at this FINAL text-column width in mm so LaTeX '
                        '\\includegraphics[width=\\textwidth] scales it ~1:1 and fonts render '
                        'at their true pt size (e.g. ~180 for an A4 1.2 cm-margin column)')
    p.add_argument('--figheight-mm', type=float, default=None, help='figure height in mm')
    p.add_argument('--point-size', type=float, default=4.0,
                   help='marker size for the panel-A UMAP scatter (was 2; raise to enlarge dots)')
    p.add_argument('--legend-marker-size', type=float, default=9.0,
                   help='size of the coloured circles in the bottom cluster legend (was 6)')
    args = p.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        w = (args.figwidth_mm / 25.4) * w_mult if args.figwidth_mm else default_w_in
        if args.figheight_mm:
            h = (args.figheight_mm / 25.4) * h_mult
        elif args.figwidth_mm:
            h = w * (default_h_in / default_w_in)  # width-only: preserve authored aspect (avoid skew)
        else:
            h = default_h_in
        return (w, h)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    print(f'reading {args.h5ad}')
    adata = sc.read_h5ad(args.h5ad)
    obs = adata.obs

    # ---- inputs / sanity checks ----
    for key in (args.cluster_key, args.pt_key, args.timepoint_key):
        if key not in obs:
            raise SystemExit(f'required obs key missing: {key}')
    if args.umap_key not in adata.obsm:
        raise SystemExit(f'required obsm key missing: {args.umap_key}')

    has_latent = args.latent_time_key in obs
    if not has_latent:
        print(f'[info] {args.latent_time_key} not in obs — skipping panel C')

    # panel D (Jonckheere-Terpstra monotone-trend) needs the cluster order file
    cluster_order = None
    if Path(args.cluster_order).exists():
        cluster_order = load_cluster_order(args.cluster_order)
    else:
        print(f'[info] {args.cluster_order} not found — skipping panel D '
              f'(Jonckheere-Terpstra trend)')
    has_jt = cluster_order is not None

    obs[args.cluster_key] = obs[args.cluster_key].astype(str)
    # Prefer a pre-annotated numeric week column (cleanest, no parsing).
    # Fall back to digit-extraction on the timepoint string only if the
    # explicit column is absent.
    if args.week_num_key in obs:
        obs['_week'] = pd.to_numeric(obs[args.week_num_key], errors='coerce')
        print(f'[ok] using obs[{args.week_num_key!r}] as numeric week '
              f'(no parsing required)')
    else:
        obs['_week'] = obs[args.timepoint_key].map(parse_week)
        print(f'[info] obs[{args.week_num_key!r}] not found — '
              f'fell back to parsing weeks from obs[{args.timepoint_key!r}]. '
              f'Run scripts/annotate_timepoint.py first for a clean week_num '
              f'column if your timepoint strings contain non-week digits '
              f'(e.g. a strain ID like "g11").')

    rng = np.random.default_rng(0)
    idx = subsample(adata, args.subsample, rng)
    print(f'plotting {len(idx):,} of {adata.n_obs:,} cells per scatter panel')

    # ---- cluster palette (alphabetical/numeric ordered) ----
    clusters = sorted(obs[args.cluster_key].unique(),
                      key=lambda s: int(s) if s.isdigit() else 1e9)
    pal = cluster_palette(clusters)                    # tab20 fallback
    if args.palette and Path(args.palette).exists():   # override with the JSON
        user_pal = {str(k): v for k, v in json.load(open(args.palette)).items()}
        hit = [c for c in clusters if c in user_pal]
        pal.update({c: user_pal[c] for c in hit})
        print(f'[ok] palette from {args.palette} '
              f'({len(hit)}/{len(clusters)} clusters matched; '
              f'rest keep the tab20 default)')
    elif args.palette:
        print(f'[info] {args.palette} not found — using the tab20 default palette')

    # ---- layout ----
    # panel indices: A=0, B=1, then C (latent) and D (JT) if present.
    n_panels = 2 + (1 if has_latent else 0) + (1 if has_jt else 0)
    idx_C = 2 if has_latent else None
    idx_D = (2 + (1 if has_latent else 0)) if has_jt else None
    if n_panels == 4:
        # 2x2 grid reads A B / C D
        # squarer default (legend now sits at the BOTTOM, not a wide right column) so at
        # text width the block fills more page vertically instead of leaving blank margin.
        # wspace 0.18 pulls the A/C and B/D columns together (less central gap) while
        # still leaving room for B/D's y-axis labels, so the columns never overlap.
        fig, axg = plt.subplots(2, 2, figsize=figsize(9.2, 9.4),
                                gridspec_kw=dict(wspace=0.18, hspace=0.30))
        axes = axg.ravel().tolist()
    else:
        fig, axes = plt.subplots(1, n_panels,
                                  figsize=figsize(5.2 * n_panels + 0.6, 4.8),
                                  gridspec_kw=dict(wspace=0.32))
        axes = np.atleast_1d(axes).ravel().tolist()

    # ============ Panel A: UMAP coloured by pseudotime ============
    ax = axes[0]
    umap = adata.obsm[args.umap_key]
    pt_vals = obs[args.pt_key].values
    order = np.argsort(pt_vals)  # plot high-pt on top so tip is visible
    sc_a = ax.scatter(umap[order, 0], umap[order, 1],
                      c=pt_vals[order], s=args.point_size, cmap='viridis',
                      alpha=0.7, lw=0, rasterized=True)
    cbar = plt.colorbar(sc_a, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label('velocity_pseudotime', fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    ax.set_title('A — pseudotime on UMAP', fontsize=9,
                 fontweight='bold', loc='left')
    ax.set_xlabel('UMAP-1', fontsize=7); ax.set_ylabel('UMAP-2', fontsize=7)
    ax.tick_params(labelsize=7)
    # add small cluster-id labels at cluster centroids (overlay)
    for cl in clusters:
        m = (obs[args.cluster_key] == cl).values
        if m.sum() < 50:
            continue
        cx, cy = umap[m, 0].mean(), umap[m, 1].mean()
        ax.text(cx, cy, cl, fontsize=7, ha='center', va='center',
                color='black', fontweight='bold',
                bbox=dict(facecolor='white', edgecolor='none',
                          alpha=0.6, pad=0.8))

    # ============ Panel B: pt vs real-time scatter ============
    ax = axes[1]
    sub = obs.iloc[idx]
    sub_pt = sub[args.pt_key].values
    sub_wk = sub['_week'].values.astype(float)
    sub_cl = sub[args.cluster_key].values
    # jitter y so a single week doesn't pile into one horizontal line
    y_jit = sub_wk + (rng.random(len(sub_wk)) - 0.5) * 0.9
    colors = [pal[cl] for cl in sub_cl]
    ax.scatter(sub_pt, y_jit, c=colors, s=3, alpha=0.5, lw=0,
               rasterized=True)
    ax.set_xlabel('velocity_pseudotime', fontsize=7)
    ax.set_ylabel('sample week', fontsize=7)
    ax.set_title('B — pseudotime ≠ real time', fontsize=9,
                 fontweight='bold', loc='left')
    ax.tick_params(labelsize=7)
    # week tick marks (only the real weeks present in the data)
    weeks = sorted(int(w) for w in pd.unique(obs['_week'].dropna()))
    ax.set_yticks(weeks)
    ax.set_yticklabels([f'wk{w}' if w > 0 else 'preCx' for w in weeks],
                       fontsize=7)
    ax.set_xlim(-0.03, 1.03)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    # bridge reference
    bridge_mask = obs[args.cluster_key] == '17'
    if bridge_mask.sum() > 0:
        bridge_pt = obs.loc[bridge_mask, args.pt_key].mean()
        ax.axvline(bridge_pt, color='#d94343', linestyle='--', lw=0.9,
                   alpha=0.7, zorder=0)
        ax.text(bridge_pt + 0.01, max(weeks) - 0.5,
                f'bridge\npt ≈ {bridge_pt:.2f}',
                fontsize=7, color='#a01010', va='top', ha='left')
    # explain the dissociation in-figure, bottom-right corner (low week + high pt is
    # empty: early-week cells all sit at low pseudotime). textwrap.fill controls the
    # box width by character count: lower `width` -> narrower (and taller) box.
    ax.text(0.99, 0.02,
            textwrap.fill('pseudotime = transcriptional dissimilarity from AR+ root; '
                          'off-diagonal cells (high week, low pt) = retained-luminal '
                          'quiescent states (e.g. cluster 0 at wk8-12, pt ≈ 0.02)',
                          width=30),
            transform=ax.transAxes, fontsize=7, style='italic',
            color='#444', ha='right', va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=2.5))

    # ============ Panel C (optional): pt vs latent_time ============
    # In addition to the scatter, computes three direction-validation
    # correlations (vpt × lt, vpt × week, lt × week) and a per-phase
    # vpt × lt ρ table — the actually-decisive numbers for the response.
    # All ρ values are computed over the full obs, not the subsampled plot.
    if has_latent:
        ax = axes[idx_C]
        sub_lt = sub[args.latent_time_key].values
        mask = ~(pd.isna(sub_pt) | pd.isna(sub_lt))
        rho_plot, _ = spearmanr(sub_pt[mask], sub_lt[mask])
        colors = [pal[cl] for cl, ok in zip(sub_cl, mask) if ok]
        ax.scatter(sub_pt[mask], sub_lt[mask], c=colors, s=3,
                   alpha=0.5, lw=0, rasterized=True)
        ax.set_xlabel('velocity_pseudotime', fontsize=7)
        ax.set_ylabel(args.latent_time_key, fontsize=7)
        sign = '+' if rho_plot >= 0 else '−'
        ax.set_title(f'C — pt vs latent_time \n'
                     f'(ρ = {sign}{abs(rho_plot):.2f})',
                     fontsize=9, fontweight='bold', loc='left')
        ax.tick_params(labelsize=7)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        # reserve empty bands top + bottom so the stats block (top-right) and the interpretation
        # cue (bottom-left) sit CLEAR of the scatter instead of on top of the points. Widen the
        # 0.55 (top) / 0.12 (bottom) fractions if the block is tall and still clips a few dots.
        _y0, _y1 = ax.get_ylim(); _yr = _y1 - _y0
        ax.set_ylim(_y0 - 0.12 * _yr, _y1 + 0.55 * _yr)

        # ---- full-data direction-validation correlations ----
        def fmt(r):
            if pd.isna(r):
                return '   n/a'
            sgn = '+' if r >= 0 else '−'
            return f'{sgn}{abs(r):.2f}'

        all_vpt = obs[args.pt_key].values
        all_lt = obs[args.latent_time_key].values
        all_wk = obs['_week'].values.astype(float)
        good = ~(pd.isna(all_vpt) | pd.isna(all_lt) | pd.isna(all_wk))
        rho_vpt_lt, _ = spearmanr(all_vpt[good], all_lt[good])
        rho_vpt_wk, _ = spearmanr(all_vpt[good], all_wk[good])
        rho_lt_wk, _ = spearmanr(all_lt[good], all_wk[good])

        # per-phase vpt × lt correlation
        cl_str = obs[args.cluster_key].astype(str).values
        phase_per_cell = np.array([PHASE_MAP.get(c, 'other') for c in cl_str])
        per_phase_rhos = []
        for group_name, phase_set in PHASE_GROUPS.items():
            m = good & np.array([p in phase_set for p in phase_per_cell])
            n = int(m.sum())
            if n < 30:
                per_phase_rhos.append((group_name, np.nan, n))
                continue
            r, _ = spearmanr(all_vpt[m], all_lt[m])
            per_phase_rhos.append((group_name, r, n))

        # stats text block — upper-left of panel, fixed-width-style alignment
        lines = [
            'Spearman ρ (whole dataset):',
            f'  velocity_pt × latent_time   {fmt(rho_vpt_lt)}',
            f'  velocity_pt × sample week   {fmt(rho_vpt_wk)}   '
            '← chronology check',
            f'  latent_time × sample week   {fmt(rho_lt_wk)}',
            '',
            'velocity_pt × latent_time, per phase:',
        ]
        for nm, r, n in per_phase_rhos:
            lines.append(f'  {nm:<28} {fmt(r)}   (n={n:,})')
        stats_text = '\n'.join(lines)
        # top-right corner (high pt + high latent_time is the large empty region;
        # the two data clouds sit at the left column and the low-lt right blob)
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
                fontsize=4, family='monospace', va='top', ha='right',
                color='#222',
                bbox=dict(facecolor='white', edgecolor='#bbb', lw=0.5,
                          alpha=0.92, pad=3))
        # interpretation cue — recognise the global-orientation-inversion +
        # within-phase-agreement pattern explicitly when present.
        valid_phase_rhos = [r for _, r, _ in per_phase_rhos if not pd.isna(r)]
        all_phases_positive = (len(valid_phase_rhos) >= 2
                                and all(r >= 0 for r in valid_phase_rhos))
        global_inversion = (not pd.isna(rho_vpt_wk) and not pd.isna(rho_lt_wk)
                             and rho_vpt_wk > 0 and rho_lt_wk < 0)
        if global_inversion and all_phases_positive:
            cue = ('global-orientation inversion of latent_time '
                   '(vpt × week +, lt × week −)\n'
                   'within-phase orderings AGREE (all per-phase ρ ≥ 0) \n'
                   '→ disagreement is global, not local\n'
                   '→ velocity_pseudotime is the chronology-aligned pt')
        elif global_inversion:
            cue = ('velocity_pt tracks real time; latent_time runs against it\n'
                   '→ velocity_pseudotime is the chronology-aligned pt')
        else:
            cue = None
        if cue:
            ax.text(0.02, 0.02, cue,
                    transform=ax.transAxes, fontsize=5, style='italic',
                    color='#a01010', ha='left', va='bottom',
                    bbox=dict(facecolor='white', edgecolor='none',
                              alpha=0.7, pad=2))

    # ============ Panel D: pseudotime vs AR->NE cluster order (JT trend) ======
    # The panel that visualises a critical claim: pseudotime rises across
    # the phenotypic adenocarcinoma->neuroendocrine order. The order is phenotypic
    # (Fig 2D), NOT a pseudotime sort, so local reversals at the terminal end are
    # expected (e.g. the DNPC-like sink cl7 sits just below the NE-mid cl14) and
    # do NOT affect the Jonckheere-Terpstra ORDERED-trend statistic (which tolerates
    # them by design). Boxes are drawn in that fixed order, so the rise is a result,
    # not a construction. Stats are on the FULL obs and match jonckheere_trend.py.
    if has_jt:
        ax = axes[idx_D]
        excl = {str(c) for c in args.exclude_clusters}
        present = [c for c in cluster_order
                   if c in set(clusters) and c not in excl]
        pt_all = obs[args.pt_key].astype(float)
        groups = [pt_all[(obs[args.cluster_key] == c).values].dropna().to_numpy()
                  for c in present]

        # boxes in AR->NE order, coloured by cluster
        bp = ax.boxplot(groups, positions=range(len(present)), widths=0.62,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', lw=1.0),
                        whiskerprops=dict(color='#888', lw=0.7),
                        capprops=dict(color='#888', lw=0.7),
                        boxprops=dict(lw=0.5, edgecolor='#555'))
        for patch, c in zip(bp['boxes'], present):
            patch.set_facecolor(pal[c]); patch.set_alpha(0.85)
        # median trend line to make the overall rise explicit
        medians = [np.median(g) if len(g) else np.nan for g in groups]
        ax.plot(range(len(present)), medians, color='#222', lw=1.1,
                marker='o', ms=2.6, zorder=6, alpha=0.85)
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels(present, fontsize=7, rotation=45, ha='right')
        ax.set_xlabel('cluster (adenocarcinoma → neuroendocrine order)', fontsize=7)
        ax.set_ylabel('pseudotime', fontsize=7)
        ax.set_title('D — pseudotime rises along AR→NE order',
                     fontsize=9, fontweight='bold', loc='left')
        ax.tick_params(labelsize=7)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

        # ---- statistics (full data; identical recipe to jonckheere_trend.py) --
        rank = {c: i for i, c in enumerate(present)}
        cl_rank = obs[args.cluster_key].map(rank).to_numpy(dtype=float)
        keep = ~np.isnan(cl_rank) & pt_all.notna().to_numpy()
        rho_ord, _ = spearmanr(pt_all.to_numpy()[keep], cl_rank[keep])
        tau_ord, _ = kendalltau(pt_all.to_numpy()[keep], cl_rank[keep])
        res_cell = jonckheere(groups)
        # timepoint-aggregated JT: each cluster contributes its per-week medians
        # (respects independence -> a believable N and Z no one can inflate)
        agg = []
        for c in present:
            sub_c = obs[obs[args.cluster_key] == c]
            med = sub_c.groupby('_week', observed=True)[args.pt_key] \
                       .median().to_numpy()
            agg.append(med[~np.isnan(med)])
        res_agg = jonckheere(agg) if all(len(g) > 0 for g in agg) else None

        # wrapped narrow so the block stays in the empty upper-left triangle
        # (the boxes rise low-left → high-right, leaving the top-left open)
        lines = [
            'Jonckheere–Terpstra',
            'ordered-trend test',
            f'  per-cell PI = {res_cell["PI"]:.2f}',
            f'    Spearman ρ = {rho_ord:.2f}',
            f'    Kendall τ = {tau_ord:.2f}',
        ]
        if res_agg is not None:
            p_agg = res_agg['P']
            p_str = f'{p_agg:.1e}' if np.isfinite(p_agg) and p_agg > 0 else '<1e-300'
            lines += [
                f'  timepoint medians (N = {res_agg["N"]}):',
                f'    PI = {res_agg["PI"]:.2f}, Z = {res_agg["Z"]:.1f}',
                f'    one-sided p = {p_str}',
            ]
        lines.append('  cell-level p omitted')
        lines.append('  (cells not independent)')
        ax.text(0.02, 0.98, '\n'.join(lines), transform=ax.transAxes,
                fontsize=4, family='monospace', va='top', ha='left',
                color='#222',
                bbox=dict(facecolor='white', edgecolor='#bbb', lw=0.5,
                          alpha=0.92, pad=3))
        print(f'[panel D] per-cell PI={res_cell["PI"]:.4f} '
              f'(ρ={rho_ord:.3f}, τ={tau_ord:.3f})'
              + (f'; aggregated N={res_agg["N"]} PI={res_agg["PI"]:.4f} '
                 f'Z={res_agg["Z"]:.3g} p={res_agg["P"]:.3e}'
                 if res_agg is not None else '; aggregated skipped'))

    # ---- shared cluster legend on the right ----
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=pal[cl], markersize=args.legend_marker_size, label=cl)
               for cl in clusters]
    fig.legend(handles=handles, loc='lower center',
               bbox_to_anchor=(0.5, -0.015), fontsize=7,
               title='cluster', title_fontsize=7, frameon=False,
               ncol=min(len(clusters), 10), columnspacing=1.0,
               handletextpad=0.4, labelspacing=0.25)

    # ---- save ----
    # two-line title: narrower than the panel block, so tight-bbox no longer pads
    # the figure left/right to fit a long single-line title.
    fig.suptitle('Pseudotime quality control:\n'
                 'distinguishing transcriptional divergence from real time',
                 fontsize=11, fontweight='bold', y=1.03)
    pdf = outdir / f'{args.stem}.pdf'
    png = outdir / f'{args.stem}.png'
    fig.savefig(pdf, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(png, bbox_inches='tight', pad_inches=0.02, dpi=200)
    print(f'wrote {pdf}\nwrote {png}')


if __name__ == '__main__':
    main()
