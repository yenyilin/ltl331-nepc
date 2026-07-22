#!/usr/bin/env python3
"""
plot_tf_activity_heatmap.py — Fig 5B: render a TF × cluster activity matrix (decoupleR
output, already z-scored) with columns ordered along the adenocarcinoma→NE axis.

Fixes the two issues we found: the decoupleR TSV columns come out STRING-sorted
(0,1,10,11,…) which is uninterpretable; this reorders them by data/cluster_order.json so
the NE-driver block reads diagonally (PRAD top-left → NE bottom-right). Works on the
existing TSVs (canonical_nepc_tfs.tsv or tf_activity_per_cluster.tsv) — no re-run needed.

Usage:
  python plot_tf_activity_heatmap.py \
      --activity Z/decoupler_tf/canonical_nepc_tfs.tsv \
      --cluster-order data/cluster_order.json \
      --out figures --stem fig5B_decoupler --plot-format pdf png
  # restrict + custom row order:
  python plot_tf_activity_heatmap.py --activity ...tf_activity_per_cluster.tsv \
      --tfs ASCL1 NEUROD1 POU3F2 INSM1 NKX2-1 FOXA2 MSX1 AR FOXA1 NKX3-1 ...
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activity", required=True, help="TF x cluster activity TSV (index=TF)")
    ap.add_argument("--cluster-order", default="data/cluster_order.json")
    ap.add_argument("--cluster-colors", default="data/cluster_colors_18.json",
                    help="cluster->hex json; colours the x-axis cluster IDs (same as the "
                         "enrichment/proliferation/ridge panels)")
    ap.add_argument("--tfs", nargs="*", default=None, help="restrict to these TF rows (else all)")
    ap.add_argument("--order-rows-by-peak", action="store_true", default=True,
                    help="sort TF rows by the axis position of their peak cluster (diagonal)")
    ap.add_argument("--no-order-rows", dest="order_rows_by_peak", action="store_false")
    ap.add_argument("--vlim", type=float, default=2.5)
    ap.add_argument("--cmap", default="RdBu_r")
    ap.add_argument("--annotate", action="store_true", help="write the z-value in each cell")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--stem", default="fig5B_decoupler")
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

    df = pd.read_csv(args.activity, sep="\t", index_col=0)
    df.columns = df.columns.astype(str)

    # ---- column order = cluster_order axis (present cols), unknowns appended ----
    if os.path.exists(args.cluster_order):
        order = [str(c) for c in json.load(open(args.cluster_order))["order"]]
        cols = [c for c in order if c in df.columns] + \
               [c for c in df.columns if c not in order]
    else:
        cols = sorted(df.columns, key=lambda s: (0, int(s)) if s.isdigit() else (1, s))
    df = df[cols]
    pos = {c: i for i, c in enumerate(cols)}

    colors = {}
    if args.cluster_colors and os.path.exists(args.cluster_colors):
        colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()}

    # ---- rows: restrict + (optional) diagonal ordering by peak-cluster position ----
    if args.tfs:
        present = [t for t in args.tfs if t in df.index]
        missing = [t for t in args.tfs if t not in df.index]
        if missing:
            print(f"[warn] TFs not in activity matrix (score by expression instead): {missing}")
        df = df.loc[present]
    if args.order_rows_by_peak and len(df):
        peak = df.to_numpy().argmax(axis=1)
        df = df.iloc[np.argsort(peak, kind="stable")]

    # ---- plot ----
    plt = setup_style(font_size=7)
    nr, nc = df.shape
    fig, ax = plt.subplots(figsize=figsize(max(5.0, 0.42 * nc + 2.0), max(3.0, 0.30 * nr + 1.5)))
    im = ax.imshow(df.to_numpy(), aspect="auto", cmap=args.cmap,
                   vmin=-args.vlim, vmax=args.vlim, interpolation="nearest")
    ax.set_xticks(range(nc)); ax.set_xticklabels(df.columns, fontsize=6, rotation=0)
    if colors:                                  # colour the cluster IDs, like the other panels
        for t, c in zip(ax.get_xticklabels(), df.columns):
            t.set_color(colors.get(str(c), "black"))
    ax.set_yticks(range(nr)); ax.set_yticklabels(df.index, fontsize=6.5)
    ax.set_xlabel("cluster (adenocarcinoma → neuroendocrine)")
    ax.set_ylabel("transcription factor")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    if args.annotate:
        for i in range(nr):
            for j in range(nc):
                v = df.iat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=4.5,
                            color="white" if abs(v) >= args.vlim * 0.6 else "#333")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("TF activity (z across clusters)", fontsize=8); cb.ax.tick_params(labelsize=6)
    ax.set_title("decoupleR TF activity along the adeno→NE axis", fontsize=9)

    os.makedirs(args.out, exist_ok=True)
    df.to_csv(os.path.join(args.out, f"{args.stem}_ordered.tsv"), sep="\t")
    savefig(fig, args.out, args.stem, args.plot_format)
    print(f"[ok] {nr} TFs x {nc} clusters; columns ordered adeno→NE; "
          f"wrote {args.stem}.{{{','.join(args.plot_format)}}} (+ _ordered.tsv) to {args.out}/")


if __name__ == "__main__":
    main()
