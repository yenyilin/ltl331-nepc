#!/usr/bin/env python3
"""
plot_scenic_decoupler_consensus.py — Supplementary Fig. S14.

Renders the SCENIC <-> decoupleR cross-method consensus (from
scenic_decoupler_consensus.py -> consensus_canonical_tfs.tsv) as a canonical-TF x
cluster map: cell colour = decoupleR TF-activity z-score, overlaid symbols mark the
per-(TF,cluster) agreement status.
The canonical NEPC drivers SCENIC did not flag as cluster-specific regulons are
independently recovered by full-cell decoupleR (decoupleR-only), while shared
calls (e.g. ASCL1 @ cluster 10) are CONCORDANT so the nominations are
method-robust, not a SCENIC subsampling artefact.

Input (wide): TF, then per cluster {clN_scenic_rank, clN_decoupler_z, clN_status}.

Usage
  python plot_scenic_decoupler_consensus.py \
      --consensus consensus_tf/consensus_canonical_tfs.tsv \
      --cluster-order data/cluster_order.json \
      --out figures/supp/S14_scenic_decoupler_consensus --plot-format pdf png
"""
import argparse
import json
import os
import re
import sys
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _utils import setup_style
    setup_style()
except Exception:
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                                "image.interpolation": "nearest"})

# status -> (marker, facecolor, edgecolor, size, label)
STATUS_STYLE = {
    "CONCORDANT":     ("o", "black",   "black", 78, "CONCORDANT (both methods)"),
    "decoupleR-only": ("D", "none",    "black", 50, "decoupleR-only (recovers what SCENIC missed)"),
    "SCENIC-only":    ("^", "none",    "0.35",  56, "SCENIC-only"),
    "DIVERGENT":      ("X", "#E45756", "black", 78, "DIVERGENT (directions disagree)"),
}


def load_cluster_order(path):
    """Return an ordered list of cluster-id strings, or None."""
    if not path or not os.path.isfile(path):
        return None
    try:
        obj = json.load(open(path))
    except Exception:
        return None
    seq = None
    if isinstance(obj, list):
        seq = obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, list) and v:
                seq = v
                break
    if not seq:
        return None
    out = []
    for x in seq:
        m = re.search(r"\d+", str(x))
        if m:
            out.append(m.group())
    return out or None


def parse_consensus(path):
    df = pd.read_csv(path, sep="\t")
    df = df.set_index("TF")
    clusters = sorted({m.group(1)
                       for c in df.columns
                       for m in [re.match(r"cl(\d+)_scenic_rank$", c)] if m},
                      key=int)
    z = pd.DataFrame(index=df.index, columns=clusters, dtype=float)
    status = pd.DataFrame(index=df.index, columns=clusters, dtype=object)
    rank = pd.DataFrame(index=df.index, columns=clusters, dtype=float)
    for c in clusters:
        z[c] = pd.to_numeric(df.get(f"cl{c}_decoupler_z"), errors="coerce")
        rank[c] = pd.to_numeric(df.get(f"cl{c}_scenic_rank"), errors="coerce")
        s = df.get(f"cl{c}_status")
        status[c] = s.astype(str).str.strip() if s is not None else ""
    return z, status, rank, clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus", default="consensus_tf/consensus_canonical_tfs.tsv")
    ap.add_argument("--cluster-order", default="data/cluster_order.json")
    ap.add_argument("--out", default="figures/supp/S14_scenic_decoupler_consensus")
    ap.add_argument("--show-rank", action="store_true",
                    help="print SCENIC rank in CONCORDANT/SCENIC-only cells")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    print(f"[load] {args.consensus}")
    z, status, rank, clusters = parse_consensus(args.consensus)

    order = load_cluster_order(args.cluster_order)
    if order:
        cols = [c for c in order if c in clusters] + \
               [c for c in clusters if c not in order]
    else:
        cols = clusters
    z, status, rank = z[cols], status[cols], rank[cols]
    tfs = list(z.index)
    print(f"[info] {len(tfs)} TFs x {len(cols)} clusters; "
          f"cluster order = {'canonical' if order else 'numeric'}")

    Z = z.values.astype(float)
    vmax = np.nanpercentile(np.abs(Z), 98)
    vmax = float(vmax) if np.isfinite(vmax) and vmax > 0 else 3.0

    fig_w = max(7.0, 0.62 * len(cols) + 2.0)
    fig_h = max(5.0, 0.42 * len(tfs) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(np.ma.masked_invalid(Z), cmap="RdBu_r",
                   norm=Normalize(vmin=-vmax, vmax=vmax),
                   aspect="auto", interpolation="nearest")
    im.cmap.set_bad("0.92")

    # gridlines between cells
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(tfs), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.0)
    ax.tick_params(which="minor", length=0)

    # status symbols
    counts = {}
    for i, tf in enumerate(tfs):
        for j, c in enumerate(cols):
            st = status.iloc[i, j]
            if st not in STATUS_STYLE:
                continue
            counts[st] = counts.get(st, 0) + 1
            mk, fc, ec, sz, _ = STATUS_STYLE[st]
            ax.scatter(j, i, marker=mk, s=sz,
                       facecolors=fc, edgecolors=ec, linewidths=1.2, zorder=3)
            if args.show_rank and st in ("CONCORDANT", "SCENIC-only"):
                rv = rank.iloc[i, j]
                if np.isfinite(rv):
                    ax.text(j + 0.34, i - 0.30, f"{int(rv)}", fontsize=9,
                            color="0.25", ha="left", va="top", zorder=4)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([str(c) for c in cols], fontsize=14, rotation=0)
    ax.set_yticks(range(len(tfs)))
    ax.set_yticklabels(tfs, fontsize=14, style="italic")
    ax.set_xlabel("cluster (adenocarcinoma → neuroendocrine order)", fontsize=15)
    ax.set_title("SCENIC ↔ decoupleR consensus for canonical NEPC regulators",
                 fontsize=17, fontweight="bold")
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cb.set_label("decoupleR TF activity (z)", fontsize=14)
    cb.ax.tick_params(labelsize=13)

    handles = [Line2D([0], [0], marker=mk, linestyle="none", markerfacecolor=fc,
                      markeredgecolor=ec, markersize=np.sqrt(sz),
                      label=textwrap.fill(f"{lab} (n={counts.get(st, 0)})", width=22))
               for st, (mk, fc, ec, sz, lab) in STATUS_STYLE.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=13, frameon=False, title="agreement", title_fontsize=14)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for fmt in args.plot_format:
        path = f"{args.out}.{fmt}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"[ok] wrote {path}")
    plt.close(fig)
    print(f"[summary] status counts: {counts}")


if __name__ == "__main__":
    main()
