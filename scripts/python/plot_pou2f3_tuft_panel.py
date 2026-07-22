#!/usr/bin/env python3
"""
Supplementary panel: the LTL331 neuroendocrine bifurcation is ASCL1+ vs a
tuft-NEGATIVE ASCL1- fate — distinct from the ASCL1-vs-tuft (POU2F3/ASCL2) split
reported in Brady 2021 (mouse GEMM), Chen 2023 (PARCB), and Ireland 2025 (SCLC).

Dot plot across all clusters (canonical adeno->NE order):
  rows  = ASCL1 driver | NE identity (CHGA/SYP/INSM1) | tuft program (POU2F3/ASCL2/OVOL3/TRPM5/AVIL)
  size  = % cells expressing
  color = mean log-normalized expression, row-scaled (per-gene min-max) for visibility
The ASCL1+ fate (cluster 10) and ASCL1- fate (clusters 1,7,9; 7=terminal) are highlighted.

Reads data/pou2f3_bifurcation_percluster.tsv (from check_pou2f3_bifurcation.py) — no h5ad
needed, runs anywhere. cluster order + stage bands from data/cluster_order.json.

Usage:
    python scripts/plot_pou2f3_tuft_panel.py
    python scripts/plot_pou2f3_tuft_panel.py --tsv data/pou2f3_bifurcation_percluster.tsv \
        --out figures/figS_pou2f3_tuft
"""
import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# gene rows, grouped (label -> genes)
GROUPS = [
    ("NE driver",      ["ASCL1"]),
    ("NE identity",    ["CHGA", "SYP", "INSM1"]),
    ("Tuft program\n(2nd fate in\nBrady/Chen/Ireland)", ["POU2F3", "ASCL2", "OVOL3", "TRPM5", "AVIL"]),
]
ASCL1POS = {"10"}
ASCL1NEG = {"1", "7", "9"}
TERMINAL = "7"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="data/pou2f3_bifurcation_percluster.tsv")
    ap.add_argument("--order", default="data/cluster_order.json")
    ap.add_argument("--out", default="figures/figS_pou2f3_tuft")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t", dtype={"cluster": str}).set_index("cluster")

    # column order
    order = None
    if os.path.exists(args.order):
        od = json.load(open(args.order))
        order = [c for c in od.get("order", []) if c in df.index]
        stages = od.get("stages", {})
    if not order:
        order = sorted(df.index, key=lambda x: (len(x), x)); stages = {}
    df = df.loc[order]

    genes = [g for _, gs in GROUPS for g in gs if f"{g}_mean" in df.columns]
    ngene, ncl = len(genes), len(order)

    # matrices: mean (row-scaled) for color, pct for size
    mean = np.array([[df.loc[c, f"{g}_mean"] for c in order] for g in genes], float)
    pct  = np.array([[df.loc[c, f"{g}_pct"]  for c in order] for g in genes], float)
    rowscaled = np.zeros_like(mean)
    for i in range(ngene):
        lo, hi = mean[i].min(), mean[i].max()
        rowscaled[i] = (mean[i] - lo) / (hi - lo) if hi > lo else 0.0

    fig, ax = plt.subplots(figsize=(0.42 * ncl + 3.2, 0.52 * ngene + 2.4))
    cmap = plt.get_cmap("Reds")
    smax = max(pct.max(), 1.0)
    for i in range(ngene):          # rows top-to-bottom
        y = ngene - 1 - i
        for j in range(ncl):
            size = 20 + (pct[i, j] / smax) * 380
            ax.scatter(j, y, s=size, color=cmap(0.15 + 0.85 * rowscaled[i, j]),
                       edgecolor="0.35", linewidth=0.4, zorder=3)

    # axes cosmetics
    ax.set_xlim(-0.8, ncl - 0.2); ax.set_ylim(-0.8, ngene - 0.2)
    ax.set_xticks(range(ncl))
    ax.set_yticks(range(ngene)); ax.set_yticklabels(genes[::-1], fontstyle="italic", fontsize=10)
    # x tick labels: highlight the fate clusters
    xtl = ax.set_xticklabels(order, fontsize=9)
    for t, c in zip(xtl, order):
        if c in ASCL1POS: t.set_color("#1f77b4"); t.set_fontweight("bold")
        elif c in ASCL1NEG: t.set_color("#d62728"); t.set_fontweight("bold")
    ax.set_xlabel("cluster (adenocarcinoma → neuroendocrine axis)", fontsize=10)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)

    # group separator lines between gene groups
    boundaries = np.cumsum([len(gs) for _, gs in GROUPS if any(f"{g}_mean" in df.columns for g in gs)])
    for b in boundaries[:-1]:
        ax.axhline(ngene - b - 0.5, color="0.85", lw=0.8, zorder=1)
    # group labels on the right
    start = 0
    for label, gs in GROUPS:
        gs = [g for g in gs if f"{g}_mean" in df.columns]
        if not gs: continue
        mid = ngene - (start + len(gs) / 2) + 0.0 - 0.5
        ax.text(ncl - 0.3, mid, label, fontsize=8, va="center", ha="left", color="0.3")
        start += len(gs)

    # highlight fate columns with faint vertical bands
    for c in list(ASCL1POS) + list(ASCL1NEG):
        j = order.index(c)
        ax.axvspan(j - 0.45, j + 0.45,
                   color=("#1f77b4" if c in ASCL1POS else "#d62728"), alpha=0.06, zorder=0)

    # title
    ax.set_title("Neuroendocrine bifurcation in LTL331: ASCL1+ vs tuft-negative ASCL1-\n"
                 "tuft markers (POU2F3/ASCL2/OVOL3) are <1% in the ASCL1- compartment",
                 fontsize=11, pad=12)

    # legends: size + color + fate colors
    size_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor="0.5",
                    markeredgecolor="0.35", markersize=np.sqrt(20 + f/smax*380)/2.5,
                    label=f"{int(f)}%") for f in [1, 25, 50, 75]]
    leg1 = ax.legend(handles=size_handles, title="% expressing", loc="upper left",
                     bbox_to_anchor=(1.02, 1.0), fontsize=8, title_fontsize=8, frameon=False,
                     labelspacing=1.1, borderpad=0.8)
    ax.add_artist(leg1)
    fate_handles = [
        Line2D([0],[0], marker="s", color="w", markerfacecolor="#1f77b4", markersize=9,
               label="ASCL1+ fate (cl 10)"),
        Line2D([0],[0], marker="s", color="w", markerfacecolor="#d62728", markersize=9,
               label="ASCL1- fate (cl 1,7,9; 7=terminal)"),
    ]
    ax.legend(handles=fate_handles, loc="lower left", bbox_to_anchor=(1.02, 0.30),
              fontsize=8, frameon=False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cax = fig.add_axes([0.995, 0.12, 0.012, 0.16])
    cb = fig.colorbar(sm, cax=cax); cb.set_label("mean expr (row-scaled)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.tight_layout(rect=[0, 0, 0.86, 1])
    for ext in ("pdf", "png"):
        fig.savefig(f"{args.out}.{ext}", dpi=200, bbox_inches="tight")
    print(f"wrote {args.out}.pdf / .png")

if __name__ == "__main__":
    main()
