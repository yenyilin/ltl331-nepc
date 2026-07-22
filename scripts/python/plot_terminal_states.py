#!/usr/bin/env python3
"""
plot_terminal_states.py — Fig 4D: the CellRank2 terminal states on the UMAP.

The point of the panel: the NE compartment has TWO terminals — ASCL1+ (cluster 10)
and ASCL1- (clusters 1/7/9, cl7 the canonical sink) — not one. CellRank's automatic
predict_terminal_states OVER-SELECTS on this discrete cluster graph: obs['term_states_fwd']
= {0, 1, 4_1, 4_2, 6, 7, 9, 10}, i.e. it also flags upstream/non-NE macrostates {0,4_1,4_2,6}
as "terminal". This script renders the curated view ( the two NE terminals highlighted )
and (with --show-dropped) draws the over-selected non-NE macrostates as faint open circles
so the curation is transparent rather than hidden.

Colours match plot_fate_umap.py (Fig 4E/F): ASCL1+ = #e7298a, ASCL1- = #1b6969, so the
terminal definition reads identically across the whole Fig 4 fate block.

Two modes:
  --by macrostate (default)  highlight the CellRank terminal-MACROSTATE cells from
                             obs['term_states_fwd'] (~30 committed cells per terminal).
                             This is the literal "terminal states" and supports --show-dropped.
  --by cluster               highlight the whole ASCL1+/- terminal CLUSTERS (fallback when
                             the object has no term_states_fwd; same view as
                             plot_fate_umap.py --terminals-panel).

Usage:
  python plot_terminal_states.py --h5ad velocity.h5ad \\
      --umap-key X_umap --cluster-key seurat_clusters --term-key term_states_fwd \\
      --ascl1-pos 10 --ascl1-neg 1 7 9 --show-dropped \\
      --colors data/cluster_colors_18.json --out figures --stem fig4D_terminal_states
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402

POS_COL, NEG_COL = "#e7298a", "#1b6969"   # ASCL1+ / ASCL1- — identical to plot_fate_umap.py
DROP_COL, BG_COL = "#b8b8b8", "#e8e8e8"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True, help="velocity / CellRank object (UMAP + term key)")
    ap.add_argument("--umap-key", default="X_umap", help="obsm key (or 'velocity_umap')")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--term-key", default="term_states_fwd",
                    help="obs column of terminal-macrostate labels (NaN for non-terminal cells)")
    ap.add_argument("--by", choices=["macrostate", "cluster"], default="macrostate",
                    help="macrostate: highlight term_states_fwd cells; cluster: whole "
                         "ASCL1+/- clusters (fallback when no terminal key)")
    ap.add_argument("--ascl1-pos", nargs="+", default=["10"],
                    help="cluster id(s) of the ASCL1+ terminal (default 10)")
    ap.add_argument("--ascl1-neg", nargs="+", default=["1", "7", "9"],
                    help="cluster id(s) of the ASCL1- terminal (default 1 7 9; matches "
                         "plot_fate_umap.py / the absorption analysis)")
    ap.add_argument("--show-dropped", action="store_true",
                    help="(--by macrostate) draw the over-selected non-NE terminal macrostates "
                         "as faint open circles, to show what predict_terminal_states picked up")
    ap.add_argument("--by-cluster-color", action="store_true",
                    help="colour each terminal by its --colors palette entry instead of the "
                         "two-group ASCL1+/- colours")
    ap.add_argument("--colors", default="data/cluster_colors_18.json",
                    help="flat {cluster: hex} palette (used by --by-cluster-color + labels)")
    ap.add_argument("--no-annotate", dest="annotate", action="store_false",
                    help="do not label each terminal with its cluster id at its centroid")
    ap.add_argument("--label-offset", nargs=2, type=float, default=[0.0, 11.0],
                    metavar=("DX", "DY"),
                    help="shift each cluster-id label by (dx, dy) in POINTS off its centroid so it "
                         "does not sit on the terminal-state dots (default 0 11 = up; try 8 8)")
    ap.add_argument("--point-size", type=float, default=7.0, help="terminal-cell marker size")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--stem", default="fig4D_terminal_states")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    ap.add_argument("--figwidth-mm", type=float, default=None,
                    help="author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa")
    ap.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    return ap.parse_args()


def main():
    args = parse_args()

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        w = (args.figwidth_mm / 25.4) * w_mult if args.figwidth_mm else default_w_in
        if args.figheight_mm:
            h = (args.figheight_mm / 25.4) * h_mult
        elif args.figwidth_mm:
            h = w * (default_h_in / default_w_in)  # width-only: preserve authored aspect (avoid skew)
        else:
            h = default_h_in
        return (w, h)

    import anndata as ad
    plt = setup_style(font_size=12)

    adata = ad.read_h5ad(args.h5ad)
    if args.umap_key not in adata.obsm:
        raise SystemExit(f"[abort] {args.umap_key} not in obsm ({list(adata.obsm)})")
    um = np.asarray(adata.obsm[args.umap_key])[:, :2]
    cl = adata.obs[args.cluster_key].astype(str).to_numpy()

    palette = {}
    if args.colors and os.path.exists(args.colors):
        palette = {str(k): v for k, v in json.load(open(args.colors)).items()}
    elif args.by_cluster_color or args.annotate:
        print(f"[warn] --colors '{args.colors}' not found; falling back to group colours")

    pos_set = [str(c) for c in args.ascl1_pos]
    neg_set = [str(c) for c in args.ascl1_neg]

    # base[i] = the cluster a cell is assigned to as a terminal; valid[i] = it IS terminal.
    if args.by == "macrostate":
        if args.term_key not in adata.obs.columns:
            raise SystemExit(
                f"[abort] obs['{args.term_key}'] not found — this object has no CellRank "
                f"terminal macrostates. Re-run with --by cluster to highlight whole "
                f"ASCL1+/- clusters instead. (obs has: {list(adata.obs.columns)[:12]} ...)")
        term = adata.obs[args.term_key]
        valid = term.notna().to_numpy()
        # macrostate label -> base cluster (e.g. '4_1' -> '4', '10' -> '10')
        base = term.astype(str).map(lambda s: s.split("_")[0]).to_numpy()
        base = np.where(valid, base, "")
    else:  # cluster mode: every cell of the pos/neg clusters is "terminal"
        base = cl.copy()
        valid = np.isin(cl, pos_set + neg_set)

    is_pos = valid & np.isin(base, pos_set)
    is_neg = valid & np.isin(base, neg_set)
    is_drop = valid & ~(is_pos | is_neg)        # over-selected non-NE terminals (macrostate mode)
    bg = ~(is_pos | is_neg | (is_drop if args.show_dropped else False))

    def colors_for(mask):
        if args.by_cluster_color and palette:
            return [palette.get(b, "#777777") for b in base[mask]]
        return None

    # ---- plot ----
    fig, ax = plt.subplots(figsize=figsize(6.0, 5.4))
    ax.scatter(um[bg, 0], um[bg, 1], s=args.point_size * 0.5, c=BG_COL, lw=0,
               rasterized=True, zorder=1)
    if args.show_dropped and is_drop.any():
        ax.scatter(um[is_drop, 0], um[is_drop, 1], s=args.point_size, facecolors="none",
                   edgecolors=DROP_COL, linewidths=0.5, rasterized=True, zorder=2,
                   label=f"over-selected (dropped): {{{','.join(sorted(set(base[is_drop])))}}}")
    neg_c = colors_for(is_neg) if args.by_cluster_color else NEG_COL
    pos_c = colors_for(is_pos) if args.by_cluster_color else POS_COL
    ax.scatter(um[is_neg, 0], um[is_neg, 1], s=args.point_size, c=neg_c, lw=0,
               rasterized=True, zorder=3, label=f"ASCL1-terminal {{{','.join(neg_set)}}}")
    ax.scatter(um[is_pos, 0], um[is_pos, 1], s=args.point_size, c=pos_c, lw=0,
               rasterized=True, zorder=4, label=f"ASCL1+terminal {{{','.join(pos_set)}}}")

    # centroid labels per terminal cluster
    if args.annotate:
        import matplotlib.patheffects as pe
        for b in sorted(set(base[is_pos | is_neg])):
            mb = (is_pos | is_neg) & (base == b)
            if mb.sum() == 0:
                continue
            cx, cy = um[mb, 0].mean(), um[mb, 1].mean()
            ax.annotate(f"cl{b}", (cx, cy), textcoords="offset points",
                        xytext=tuple(args.label_offset), ha="center", va="center",
                        fontsize=10, fontweight="bold", color="white", zorder=6,
                        path_effects=[pe.withStroke(linewidth=1.6, foreground="black")])

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=11); ax.set_ylabel("UMAP2", fontsize=8)
    title = ("NE terminal states\n(CellRank macrostates)" if args.by == "macrostate"
             else "NE terminal clusters")
    ax.set_title(title, fontsize=12)
    ax.legend(loc="best", fontsize=8, frameon=False, markerscale=2.0)
    for s in ax.spines.values():
        s.set_visible(False)

    os.makedirs(args.out, exist_ok=True)
    savefig(fig, args.out, args.stem, args.plot_format)
    plt.close(fig)

    # ---- summary tsv: per terminal cluster, n highlighted cells + group ----
    rows = []
    for grp, mask in (("ASCL1+", is_pos), ("ASCL1-", is_neg),
                      ("dropped", is_drop if args.show_dropped else np.zeros_like(is_pos))):
        for b in sorted(set(base[mask])):
            rows.append({"group": grp, "cluster": b, "n_cells": int((mask & (base == b)).sum())})
    pd.DataFrame(rows).to_csv(os.path.join(args.out, f"{args.stem}_counts.tsv"),
                             sep="\t", index=False)
    print(f"[ok] {args.by} mode: ASCL1+ {int(is_pos.sum())} cells / "
          f"ASCL1- {int(is_neg.sum())} cells"
          + (f" / dropped {int(is_drop.sum())}" if args.show_dropped else ""))
    print(f"[ok] wrote {args.stem}.{{{','.join(args.plot_format)}}} (+ _counts.tsv) to {args.out}/")


if __name__ == "__main__":
    main()
