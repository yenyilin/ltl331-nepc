#!/usr/bin/env python3
"""
plot_chea3_map.py — Supplementary Fig. S15 panel C: ChEA3 binding-evidence map.

A TF (rows) x cluster-signature (cols) heatmap of ChEA3 Integrated--meanRank binding rank,
the binding-evidence twin of the SCENIC (Fig 5A) and decoupleR (Fig 5B) regulon maps.
Dark = top binding-supported (low rank). The regulon-only TFs (MSX1/PAX7/RXRG) are shown as a
greyed block to make explicit that they are NOT binding-supported (they are co-expression/regulon
nominations, consistent with the manuscript's "regulon activity, not expression" framing).

Reads data/chea3_input/cl<N>_chea3.json (produced by the ChEA3 API run).
Outputs figures/fig5_S15C_chea3_map.{pdf,png}.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse

_ap = argparse.ArgumentParser()
_ap.add_argument("--figwidth-mm", type=float, default=None,
                 help="author panel at this final width in mm (S15 cell target) so it scales ~1:1")
_ap.add_argument("--figheight-mm", type=float, default=None)
_args, _ = _ap.parse_known_args()

def _figsize(dw, dh):
    fw, fh = _args.figwidth_mm, _args.figheight_mm
    if fw and fh: return (fw / 25.4, fh / 25.4)
    if fw:        return (fw / 25.4, dh * (fw / 25.4) / dw)
    if fh:        return (dw * (fh / 25.4) / dh, fh / 25.4)
    return (dw, dh)

BASE = "data/chea3_input"
N_TOTAL = 1632  # TFs ranked per query (for log normalization)

# columns = trajectory order: priming -> delamination -> early NE -> ASCL1+ / ASCL1-
# (stage labels match manuscript vocabulary: cl4 "priming", cl17 "delamination",
#  cl13 "early-neuroendocrine" — NOT "specification", which would invert the neural-crest
#  tier order the paper invokes, results.md:22/36)
COLS = [("cl4",  "cl4",     "priming"),
        ("cl17", "cl17",    "delamination"),
        ("cl13", "cl13",    "early NE"),
        ("cl10", "cl10",    "ASCL1+ NE"),
        ("cl7",  "cl7_NE",  "ASCL1- NE")]

# row groups (label, [TFs], is_regulon_only)
GROUPS = [
    ("EMT/delam.",   ["TWIST1","TWIST2","SNAI2","ZEB1","ZEB2","PRRX1"], False),
    ("NE drivers",   ["ASCL1","NEUROD1","POU3F2","INSM1","SOX2"],       False),
    ("NE lin.",      ["MYT1L","ST18"],                                  False),
    ("regulon",      ["MSX1","PAX7","RXRG"],                            True),
]

def load_ranks(stem):
    d = json.load(open(f"{BASE}/{stem}_chea3.json"))
    return {r["TF"]: int(r["Rank"]) for r in d["Integrated--meanRank"]}

ranks = {col_key: load_ranks(stem) for col_key, stem, _ in COLS}

rows, row_groups, regulon_mask = [], [], []
for glabel, tfs, ro in GROUPS:
    for tf in tfs:
        rows.append(tf); row_groups.append(glabel); regulon_mask.append(ro)

nrow, ncol = len(rows), len(COLS)
strength = np.full((nrow, ncol), np.nan)   # log10(N/rank); dark = top-ranked
rankmat  = np.full((nrow, ncol), -1, int)
for j, (col_key, _, _) in enumerate(COLS):
    for i, tf in enumerate(rows):
        rk = ranks[col_key].get(tf)
        if rk:
            rankmat[i, j] = rk
            strength[i, j] = np.log10(N_TOTAL / rk)

ro_arr = np.array(regulon_mask)
# split: supported rows get the colormap; regulon-only rows forced grey
supp = np.where(ro_arr[:, None], np.nan, strength)

fig, ax = plt.subplots(figsize=_figsize(0.95*ncol + 3.2, 0.42*nrow + 1.6))
cmap = plt.cm.magma_r.copy(); cmap.set_bad("#f0f0f0")
vmax = np.nanmax(strength)
im = ax.imshow(np.ma.masked_invalid(supp), cmap=cmap, vmin=0, vmax=vmax, aspect="auto")

# grey band + grey cells for the regulon-only block
for i in range(nrow):
    if ro_arr[i]:
        ax.add_patch(Rectangle((-0.5, i-0.5), ncol, 1, color="#d9d9d9", zorder=1))

# annotate each cell with its rank
for i in range(nrow):
    for j in range(ncol):
        rk = rankmat[i, j]
        txt = str(rk) if rk > 0 else "—"
        if ro_arr[i]:
            color = "#444444"
        else:
            s = strength[i, j]
            color = "white" if (not np.isnan(s) and s > vmax*0.55) else "black"
        ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color, zorder=3)

# ticks
ax.set_xticks(range(ncol))
_stg = {'priming': 'prime', 'delamination': 'delam', 'early NE': 'early-NE',
        'ASCL1+ NE': 'ASCL1+', 'ASCL1- NE': 'ASCL1−'}
ax.set_xticklabels([f"{c[0]} · {_stg.get(c[2], c[2])}" for c in COLS],
                   fontsize=7, rotation=30, ha='right', rotation_mode='anchor')
ax.set_yticks(range(nrow)); ax.set_yticklabels(rows, fontsize=8)
ax.set_xlim(-0.5, ncol-0.5); ax.set_ylim(nrow-0.5, -0.5)

# row-group brackets/labels on the left
start = 0
for glabel, tfs, ro in GROUPS:
    end = start + len(tfs)
    ax.add_patch(Rectangle((-0.5, start-0.5), ncol, len(tfs), fill=False,
                           edgecolor="#888888" if not ro else "#aaaaaa", lw=1.2, zorder=2))
    ax.text(ncol-0.4, (start+end-1)/2, glabel, rotation=270, va="center", ha="left",
            fontsize=6.5, color="#555555")
    if start > 0:
        ax.axhline(start-0.5, color="white", lw=2)
    start = end

cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.16)
cbar.set_label("ChEA3 binding support\nlog₁₀(N / rank)  (dark = top-ranked)", fontsize=7.5)
ax.set_title("ChEA3 TF-target binding rank", fontsize=9)
plt.tight_layout()
os.makedirs("figures", exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(f"figures/fig5_S15C_chea3_map.{ext}", dpi=200, bbox_inches="tight")
print("wrote figures/fig5_S15C_chea3_map.{pdf,png}")
print(f"  rows={nrow} (incl {int(ro_arr.sum())} regulon-only), cols={ncol}")
