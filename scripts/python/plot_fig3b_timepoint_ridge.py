#!/usr/bin/env python3
"""
plot_fig3b_timepoint_ridge.py — Fig 3B: per-timepoint ridgeline of the cancer-cell
distribution along the adenocarcinoma->neuroendocrine phenotypic axis, scaled per
timepoint, with the four stage bands (AR-high / AR-low / Intermediate / NEPC) as the
x-axis background. Reproduces the original Seurat ridge panel in scanpy/matplotlib.

x-axis (--x-axis):
  score   (DEFAULT, recommended) — each cell's continuous AR→NE phenotype score (NE − PRAD,
          scored on .raw or from --ne-key/--prad-key). Monotonic in time and free of the
          within-stage ordering artifact: the AR-low timepoints (wk2–12) correctly stack in
          the AR-low band instead of being spread by arbitrary cluster index, and later
          timepoints can't appear "earlier". Stage bands = score ranges split at the midpoints
          of each stage's median score (asserts NO within-stage cluster ordering).
  cluster — the canonical-cluster-index axis (jittered). Kept for visual comparison; note the
          numeric within-stage order {0,2,5,8,12,16} does NOT track time, so a timepoint's peak
          can land left of an earlier one (the artifact that motivated the score axis).
One ridge per timepoint (rows), each KDE normalized to its own max (the "scaled to even out
representation" of the caption), stacked from early (bottom) to late (top).

Uses obs only (cluster + timepoint) — no expression, fast. The 7 cluster_order
stages are collapsed to the 4 caption stages (AR-high; AR-low incl. pre-EMT;
Intermediate; NEPC) for the background bands.

Usage
  python plot_fig3b_timepoint_ridge.py --h5ad ltl331_base.h5ad \\
      --cluster-key seurat_clusters --time-key week_num \\
      --cluster-order data/cluster_order.json --outdir figures --stem fig3B
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig, numkey  # noqa: E402

# collapse the 7 cluster_order stages -> the 4 caption stages, with band colours.
# Colours follow the canonical phenotype palette (data/cluster_phenotype_colors.json /
# plot_timepoint_panel.py PHENO_COLOR) so the AR-high->NEPC mapping is purple->orange->
# magenta->green across every Fig 2/3 panel — NOT the old blue/teal/orange/red here.
FOUR_STAGE = [
    ("AR-high",      ("#7570b3", 0.18), lambda s: "ar-high" in s),
    ("AR-low",       ("#d95f02", 0.18), lambda s: ("ar-low" in s) or ("pre-emt" in s)),
    ("Intermediate", ("#e7298a", 0.22), lambda s: "intermediate" in s),
    ("NEPC",         ("#1b9e77", 0.16), lambda s: s.startswith("ne") or "nepc" in s),
]


# default phenotype-score gene sets (used by --x-axis score if no obs keys given)
NE_DEFAULT = ["NCAM1", "CHGA", "CHGB", "SYP", "ENO2", "ASCL1", "INSM1"]
PRAD_DEFAULT = ["AR", "KLK3", "KLK2", "FKBP5", "TMPRSS2", "NKX3-1", "FOXA1"]


def stage_bands(stages, pos):
    """CLUSTER mode: stages {name:[ids]}, pos {cluster:x} -> [(label,color,a,x0,x1)] by index."""
    bands = []
    for label, (color, alpha), match in FOUR_STAGE:
        xs = [pos[str(c)] for sname, cls in stages.items() if match(sname.lower())
              for c in cls if str(c) in pos]
        if xs:
            bands.append((label, color, alpha, min(xs) - 0.5, max(xs) + 0.5))
    return bands


def collapse_stage(stages):
    """{cluster id -> 4-stage label} via FOUR_STAGE matchers."""
    out = {}
    for sname, cls in stages.items():
        for label, _, match in FOUR_STAGE:
            if match(sname.lower()):
                for c in cls:
                    out[str(c)] = label
                break
    return out


def score_bands(stages, cl_arr, xvals):
    """SCORE mode: contiguous bands partitioned at midpoints of each stage's median score,
    so the x-axis is split AR-high<...<NEPC by actual phenotype (no within-stage ordering)."""
    s_of = collapse_stage(stages)
    colormap = {label: (c, a) for label, (c, a), _ in FOUR_STAGE}
    meds = {}
    for label in colormap:
        m = np.array([s_of.get(str(c)) == label for c in cl_arr])
        if m.any():
            meds[label] = float(np.median(xvals[m]))
    if not meds:
        return []
    ordered = sorted(meds, key=lambda l: meds[l])          # by score, low->high
    cents = [meds[l] for l in ordered]
    xmin, xmax = float(xvals.min()), float(xvals.max())
    edges = [xmin] + [(cents[i] + cents[i + 1]) / 2 for i in range(len(cents) - 1)] + [xmax]
    return [(l, colormap[l][0], colormap[l][1], edges[i], edges[i + 1])
            for i, l in enumerate(ordered)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--time-key", default="week_num",
                    help="obs column with the timepoint (numeric weeks preferred)")
    ap.add_argument("--cluster-order", default="data/cluster_order.json")
    ap.add_argument("--cluster-colors", default="data/cluster_colors_18.json",
                    help="cluster->hex json; colours the x-axis cluster IDs in --x-axis cluster "
                         "mode (same as the enrichment/proliferation panels)")
    ap.add_argument("--exclude-clusters", nargs="*", default=["18"])
    ap.add_argument("--x-axis", choices=["score", "cluster"], default="score",
                    help="'score' (default, RECOMMENDED) = continuous AR→NE phenotype score "
                         "(NE−PRAD), monotonic in time, no within-stage ordering artifact; "
                         "'cluster' = the old canonical-cluster-index axis (for visual comparison).")
    ap.add_argument("--ne-key", default=None, help="obs column with a precomputed NE score")
    ap.add_argument("--prad-key", default=None, help="obs column with a precomputed PRAD/AR score")
    ap.add_argument("--ne-genes", nargs="+", default=NE_DEFAULT,
                    help="NE gene set if --ne-key not given (scored on .raw)")
    ap.add_argument("--prad-genes", nargs="+", default=PRAD_DEFAULT,
                    help="PRAD/AR gene set if --prad-key not given (scored on .raw)")
    ap.add_argument("--bw", type=float, default=None,
                    help="KDE bandwidth; cluster mode = factor/std (default 0.5), "
                         "score mode = scalar bw_method (default 'scott')")
    ap.add_argument("--height", type=float, default=0.9, help="ridge overlap (1=touching)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--stem", default="fig3B")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    ap.add_argument("--figwidth-mm", type=float, default=None,
                    help="author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa")
    ap.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    args = ap.parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
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

    rng = np.random.default_rng(args.seed)

    adata = sc.read_h5ad(args.h5ad)
    obs = adata.obs[[args.cluster_key, args.time_key]].copy()
    obs[args.cluster_key] = obs[args.cluster_key].astype(str)
    if args.exclude_clusters:
        obs = obs[~obs[args.cluster_key].isin([str(c) for c in args.exclude_clusters])]

    co = json.load(open(args.cluster_order))
    order = [str(c) for c in co.get("order", [])]
    present = [c for c in order if c in set(obs[args.cluster_key])]
    pos = {c: i for i, c in enumerate(present)}
    obs = obs[obs[args.cluster_key].isin(pos)]

    colors = {}
    if args.cluster_colors and os.path.exists(args.cluster_colors):
        colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()}

    if args.x_axis == "cluster":
        x = obs[args.cluster_key].map(pos).to_numpy(float) + rng.uniform(-0.45, 0.45, len(obs))
    else:  # score: continuous AR→NE phenotype score (NE − PRAD)
        cells = obs.index
        if args.ne_key and args.prad_key and args.ne_key in adata.obs and args.prad_key in adata.obs:
            ne = adata.obs.loc[cells, args.ne_key].to_numpy(float)
            pr = adata.obs.loc[cells, args.prad_key].to_numpy(float)
        else:
            src = set(map(str, (adata.raw.var_names if adata.raw is not None else adata.var_names)))
            ne_g = [g for g in args.ne_genes if g in src]
            pr_g = [g for g in args.prad_genes if g in src]
            if not ne_g or not pr_g:
                raise SystemExit(f"[abort] need NE & PRAD genes present; NE={ne_g} PRAD={pr_g}")
            print(f"[score] NE={ne_g}  PRAD={pr_g}  (scored on .raw)")
            sc.tl.score_genes(adata, ne_g, score_name="_ne_sc", use_raw=adata.raw is not None)
            sc.tl.score_genes(adata, pr_g, score_name="_prad_sc", use_raw=adata.raw is not None)
            ne = adata.obs.loc[cells, "_ne_sc"].to_numpy(float)
            pr = adata.obs.loc[cells, "_prad_sc"].to_numpy(float)
        x = ne - pr
    obs = obs.assign(_x=x)

    # timepoints sorted (numeric if possible)
    tp = obs[args.time_key]
    try:
        tkey = lambda v: float(v)
        times = sorted(tp.dropna().unique(), key=tkey)
    except (TypeError, ValueError):
        times = sorted(tp.dropna().unique(), key=lambda v: numkey(str(v)))

    def tlabel(v):
        try:
            fv = float(v)
            return "preCx" if fv == 0 else f"wk{int(fv) if fv.is_integer() else fv}"
        except (TypeError, ValueError):
            return str(v)

    from scipy.stats import gaussian_kde
    xall = obs["_x"].to_numpy()
    if args.x_axis == "cluster":
        gmin, gmax = -0.5, len(present) - 0.5
    else:
        pad = 0.04 * (xall.max() - xall.min())
        gmin, gmax = xall.min() - pad, xall.max() + pad
    grid = np.linspace(gmin, gmax, 400)
    span = gmax - gmin

    plt = setup_style(font_size=8)
    fig, ax = plt.subplots(figsize=figsize(7.5, 0.55 * len(times) + 1.8))

    # background stage bands (behind everything)
    if args.x_axis == "cluster":
        bands = stage_bands(co.get("stages", {}), pos)
    else:
        bands = score_bands(co.get("stages", {}), obs[args.cluster_key].to_numpy(), xall)
    ytop = len(times) - 1 + args.height + 0.4
    # stage labels: stagger into two vertical tiers so neighbouring labels never collide
    # when a band is narrow (e.g. the Intermediate band in score mode) and its word is wide.
    # bands run AR-high / AR-low / Intermediate / NEPC, so alternating tiers separates every
    # adjacent pair. A faint leader ties each raised label back to its band.
    stag = 0.7
    for j, (label, color, alpha, x0, x1) in enumerate(bands):
        ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0, zorder=0)
        xc, yj = (x0 + x1) / 2, ytop + (j % 2) * stag
        if j % 2:
            ax.plot([xc, xc], [ytop + 0.3, yj], color=color, lw=0.5, alpha=0.6,
                    zorder=1, clip_on=False)
        ax.text(xc, yj, label, ha="center", va="bottom",
                fontsize=6, color=color, fontweight="bold")

    # ridges, early (bottom) -> late (top)
    for i, t in enumerate(times):
        xi = obs.loc[tp == t, "_x"].to_numpy()
        if len(xi) < 5:
            continue
        if args.x_axis == "cluster":
            bw = (args.bw if args.bw is not None else 0.5) / max(xi.std(), 1e-6)
        else:
            bw = args.bw if args.bw is not None else "scott"
        d = gaussian_kde(xi, bw_method=bw)(grid)
        d = d / d.max() * args.height           # scale each ridge to its own max
        base = float(i)
        ax.fill_between(grid, base, base + d, color="0.25", alpha=0.85, lw=0, zorder=2 + i)
        ax.plot(grid, base + d, color="white", lw=0.6, zorder=2 + i)
        ax.text(gmin - 0.02 * span, base + 0.05, tlabel(t), ha="right", va="bottom", fontsize=7)

    ax.set_xlim(gmin, gmax)
    ax.set_ylim(-0.2, ytop + stag + 0.7)
    ax.set_yticks([])
    if args.x_axis == "cluster":
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels(present, fontsize=6)
        if colors:                              # colour the cluster IDs, like enrichment/proliferation
            for t, c in zip(ax.get_xticklabels(), present):
                t.set_color(colors.get(c, "black"))
        ax.set_xlabel("cluster")
    else:
        ax.axvline(0, color="0.6", lw=0.5, ls="--", zorder=1)
        ax.set_xlabel("phenotype score  \n(PRAD/AR  ←   NE−PRAD   →  NE)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    os.makedirs(args.outdir, exist_ok=True)
    savefig(fig, args.outdir, args.stem, args.plot_format)
    plt.close(fig)
    print(f"[ok] wrote {args.stem}.{{{','.join(args.plot_format)}}} to {args.outdir}/ "
          f"({len(times)} timepoints, x-axis={args.x_axis})")


if __name__ == "__main__":
    main()
