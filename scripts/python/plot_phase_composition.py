#!/usr/bin/env python3
"""
plot_phase_composition.py — cell-cycle phase (G1/S/G2M) composition as stacked bars,
per seurat_cluster and per timepoint.

For each cluster (and each timepoint) the bar shows the fraction of its cells in each
phase (sums to 1), so you can read where proliferative (S/G2M) cells concentrate along
the adenocarcinoma->NE axis and across the castration time course. Uses obs only
(phase + cluster + timepoint) — no expression matrix needed.

  --by cluster | timepoint | both (default)   which panel(s) to draw
  --counts                                     absolute cell counts instead of fractions

Clusters are ordered by data/cluster_order.json (adeno->NE); timepoints chronologically
by week_num. Phase order and colours are configurable.

Example:
  python plot_phase_composition.py --h5ad ltl331_base.h5ad \\
      --phase-key phase --cluster-key seurat_clusters --timepoint-key timepoint \\
      --cluster-order data/cluster_order.json --outdir figures --stem phase_composition
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig, numkey  # noqa: E402

# conventional proliferation ramp: G1 resting (grey) -> S (amber) -> G2M mitotic (red)
PHASE_COLORS = {"G1": "#bdbdbd", "S": "#fdae61", "G2M": "#d73027"}


def order_timepoints(obs, tp_key, week_key):
    tp = obs[tp_key].astype(str)
    if week_key in obs:
        wk = pd.to_numeric(obs[week_key], errors="coerce").groupby(tp).median()
        return list(wk.sort_values().index)
    def k(s):
        if "pre" in s.lower():
            return -1
        d = "".join(c for c in s if c.isdigit())
        return int(d) if d else 999
    return sorted(tp.unique(), key=k)


def composition(obs, row_key, row_order, phase_key, phase_order, counts):
    """rows = row_order, cols = phase_order; values = fraction (or counts)."""
    ct = pd.crosstab(obs[row_key].astype(str), obs[phase_key].astype(str))
    ct = ct.reindex(index=row_order, columns=phase_order, fill_value=0)
    if not counts:
        ct = ct.div(ct.sum(axis=1), axis=0).fillna(0)
    return ct


def draw(ax, ct, colors, counts, xlabel, title):
    x = np.arange(len(ct.index))
    bottom = np.zeros(len(ct.index))
    for ph in ct.columns:
        ax.bar(x, ct[ph].values, bottom=bottom, width=0.82,
               color=colors.get(ph, "#777777"), edgecolor="white", lw=0.3, label=ph)
        bottom += ct[ph].values
    ax.set_xticks(x); ax.set_xticklabels(ct.index, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("number of cells" if counts else "fraction of cells", fontsize=9)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_xlim(-0.6, len(ct.index) - 0.4)
    if not counts:
        ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--phase-key", default="phase")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--timepoint-key", default="timepoint")
    ap.add_argument("--week-key", default="week_num")
    ap.add_argument("--cluster-order", default="data/cluster_order.json")
    ap.add_argument("--exclude-clusters", nargs="*", default=["18"],
                    help="cluster ids to drop (default 18 = inferCNV reference)")
    ap.add_argument("--by", choices=["cluster", "timepoint", "both"], default="both")
    ap.add_argument("--layout", choices=["horizontal", "vertical"], default="horizontal",
                    help="for --by both: side-by-side (horizontal) or stacked "
                         "top=cluster / bottom=timepoint (vertical)")
    ap.add_argument("--phase-order", nargs="+", default=["G1", "S", "G2M"])
    ap.add_argument("--phase-colors", default=None, help="JSON {phase: hex} overriding defaults")
    ap.add_argument("--counts", action="store_true", help="absolute counts instead of fractions")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--stem", default="phase_composition")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    import anndata as ad
    plt = setup_style(font_size=8)
    from matplotlib.patches import Patch

    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs
    for k in (args.phase_key, args.cluster_key):
        if k not in obs.columns:
            raise SystemExit(f"[abort] '{k}' not in obs ({list(obs.columns)[:15]} ...)")
    obs = obs.copy()
    obs[args.cluster_key] = obs[args.cluster_key].astype(str)

    if args.exclude_clusters:
        drop = obs[args.cluster_key].isin([str(c) for c in args.exclude_clusters])
        if drop.any():
            print(f"excluding clusters {args.exclude_clusters}: {int(drop.sum())} cells")
            obs = obs[~drop]

    colors = dict(PHASE_COLORS)
    if args.phase_colors and os.path.exists(args.phase_colors):
        colors.update({str(k): v for k, v in json.load(open(args.phase_colors)).items()})

    # keep only phases that are present, in the requested order
    present_ph = [p for p in args.phase_order if p in set(obs[args.phase_key].astype(str))]
    extra = sorted(set(obs[args.phase_key].astype(str)) - set(args.phase_order))
    if extra:
        print(f"[warn] phases not in --phase-order (appended): {extra}")
    phase_order = present_ph + extra

    # cluster order (phenotypic) + present
    present_cl = sorted(set(obs[args.cluster_key]), key=lambda s: numkey(s))
    if args.cluster_order and os.path.exists(args.cluster_order):
        co = [str(c) for c in json.load(open(args.cluster_order)).get("order", [])]
        cl_order = [c for c in co if c in present_cl] + [c for c in present_cl if c not in co]
    else:
        cl_order = present_cl
    tp_order = order_timepoints(obs, args.timepoint_key, args.week_key) \
        if args.timepoint_key in obs.columns else None
    if args.by in ("timepoint", "both") and tp_order is None:
        raise SystemExit(f"[abort] --timepoint-key '{args.timepoint_key}' not in obs")

    panels = []
    if args.by in ("cluster", "both"):
        panels.append(("cluster", composition(obs, args.cluster_key, cl_order,
                                              args.phase_key, phase_order, args.counts),
                       "cluster (adenocarcinoma → neuroendocrine)", "Phase composition per cluster"))
    if args.by in ("timepoint", "both"):
        panels.append(("timepoint", composition(obs, args.timepoint_key, tp_order,
                                                args.phase_key, phase_order, args.counts),
                       "timepoint", "Phase composition per timepoint"))

    widths = [max(2.5, 0.32 * len(ct.index) + 1.2) for _, ct, _, _ in panels]
    handles = [Patch(facecolor=colors.get(p, "#777777"), label=p) for p in phase_order]
    if args.layout == "vertical" and len(panels) > 1:
        # stacked: top = first panel (cluster), bottom = second (timepoint)
        fig, axes = plt.subplots(len(panels), 1, figsize=(max(widths) + 1.5, 3.4 * len(panels)),
                                 squeeze=False)
        axcol = [axes[i][0] for i in range(len(panels))]
        for ax, (_, ct, xlab, title) in zip(axcol, panels):
            draw(ax, ct, colors, args.counts, xlab, title)
        axcol[0].legend(handles=handles, title="phase", loc="center left",
                        bbox_to_anchor=(1.02, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    else:
        fig, axes = plt.subplots(1, len(panels), figsize=(sum(widths) + 1.5, 4.2),
                                 gridspec_kw={"width_ratios": widths}, squeeze=False)
        for ax, (_, ct, xlab, title) in zip(axes[0], panels):
            draw(ax, ct, colors, args.counts, xlab, title)
        axes[0][-1].legend(handles=handles, title="phase", loc="center left",
                           bbox_to_anchor=(1.02, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    fig.tight_layout()

    os.makedirs(args.outdir, exist_ok=True)
    savefig(fig, args.outdir, args.stem, args.plot_format)
    plt.close(fig)
    # also write the underlying tables
    for name, ct, _, _ in panels:
        ct.to_csv(os.path.join(args.outdir, f"{args.stem}_by_{name}.tsv"), sep="\t")
    print(f"[ok] wrote {args.stem}.{{{','.join(args.plot_format)}}} (+ tables) to {args.outdir}/")


if __name__ == "__main__":
    main()
