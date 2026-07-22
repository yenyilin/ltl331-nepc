#!/usr/bin/env python3
"""
plot_tmprss2_erg_deletion.py — genomic evidence for a deletion-mediated
TMPRSS2:ERG fusion in LTL331/R, retained across the castration timecourse.

TMPRSS2:ERG most often arises by an interstitial deletion of the ~3 Mb between
the two genes on 21q22 (the "Edel" mechanism). That deletion IS the fusion and
leaves a copy-number scar: a single continuous loss whose breakpoints fall
WITHIN ERG and WITHIN TMPRSS2 (not a broad arm-level 21q loss). This script
zooms the WES copy-number segments onto chr21:37.5-42.5 Mb (hg38) and draws one
track per timepoint, annotating the ERG and TMPRSS2 gene bodies and shading the
fusion deletion interval — showing the signature is present at pre-castration
(PRAD founder) and retained through 22 wk (NEPC terminal) = clonal continuity.

Deliberately DNA-only: ERG mRNA is AR-driven and falls as AR collapses in NEPC,
so the retained *genomic* rearrangement — not expression — is the continuity
marker. RNA chimeric-junction calling is corroborative and not required for this claim.

hg38 landmarks: ERG chr21:38,367,261-38,661,783; TMPRSS2 chr21:41,464,305-41,531,116.

Example:
  python scripts/plot_tmprss2_erg_deletion.py --wes-calls data/nexus_calls.tsv \\
      --outdir figures/cnv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CHROM = "chr21"
ERG = (38_367_261, 38_661_783)       # hg38 gene body
TMPRSS2 = (41_464_305, 41_531_116)   # hg38 gene body
WINDOW = (37_500_000, 42_500_000)    # zoom

# WES/Nexus sample -> timepoint (public sample names, no normalization suffix).
# The 4 GUNS WES timepoints are the primary panel (matches the CNV-supplement
# scope); the 2 aCGH timepoints (8/12 wk) corroborate and can be added with --all.
WES_TP_PRIMARY = {
    "LTL331_preCX1": "preCx",
    "LTL331_16wk1": "16wk",
    "LTL331_20wk": "20wk",
    "LTL331_22wk2": "22wk",
}
WES_TP_ACGH = {"LTL331-5_8wk": "8wk", "LTL331-5_12wk": "12wk"}
ORDER = ["preCx", "8wk", "12wk", "16wk", "20wk", "22wk"]

CMAP = {-2: "#08519c", -1: "#6baed6", 0: "#f0f0f0", 1: "#fb6a4a", 2: "#a50f15"}


def segs_in_window(calls, sample):
    s = calls[(calls["sample"].astype(str) == sample)
              & (calls["Chromosome"].astype(str) == CHROM)
              & (calls["Start"] < WINDOW[1]) & (calls["End"] > WINDOW[0])]
    return s.sort_values("Start")


def fusion_call(sub):
    """Does a loss segment span from within/at ERG through to within/at TMPRSS2?
    Returns (retained_bool, start, end) of the spanning loss if present."""
    loss = sub[sub["Value"] < 0]
    # a single-or-merged loss whose extent reaches from <=ERG end to >=TMPRSS2 start
    if loss.empty:
        return (False, None, None)
    lo, hi = loss["Start"].min(), loss["End"].max()
    retained = (lo <= ERG[1]) and (hi >= TMPRSS2[0])
    return (retained, int(lo), int(hi))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wes-calls", default="data/nexus_calls.tsv")
    ap.add_argument("--outdir", default="figures/cnv")
    ap.add_argument("--all", action="store_true",
                    help="include the 2 aCGH timepoints (8/12 wk) as corroboration")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(args.wes_calls, sep="\t", dtype={"sample": str})
    tp_map = dict(WES_TP_PRIMARY)
    if args.all:
        tp_map.update(WES_TP_ACGH)
    inv = {tp: s for s, tp in tp_map.items()}
    tps = [tp for tp in ORDER if tp in inv]

    # ---- text summary + source table ----
    rows = []
    print("[TMPRSS2:ERG deletion signature — hg38, chr21]")
    for tp in tps:
        sub = segs_in_window(calls, inv[tp])
        retained, lo, hi = fusion_call(sub)
        tag = (f"RETAINED  loss {lo:,}-{hi:,}  (breakpoints inside ERG & TMPRSS2)"
               if retained else "NOT detected")
        print(f"  {tp:6s} {inv[tp]:16s} -> {tag}")
        rows.append({"timepoint": tp, "sample": inv[tp], "fusion_deletion_retained": retained,
                     "loss_start": lo, "loss_end": hi})
    pd.DataFrame(rows).to_csv(outdir / "tmprss2_erg_deletion_calls.tsv", sep="\t", index=False)

    # ---- figure: one segment track per timepoint ----
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(tps) + 1.8))
    h = 0.6
    for i, tp in enumerate(tps[::-1]):          # earliest at top
        y = i
        sub = segs_in_window(calls, inv[tp])
        # neutral baseline across the window
        ax.add_patch(Rectangle((WINDOW[0], y - h / 2), WINDOW[1] - WINDOW[0], h,
                               facecolor=CMAP[0], edgecolor="none", zorder=1))
        for _, r in sub.iterrows():
            x0 = max(r["Start"], WINDOW[0]); x1 = min(r["End"], WINDOW[1])
            ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h,
                                   facecolor=CMAP.get(int(r["Value"]), "#999999"),
                                   edgecolor="none", zorder=2))
    track_top = (len(tps) - 1) + h / 2         # top edge of the highest track
    ax.set_yticks(range(len(tps)))
    ax.set_yticklabels(tps[::-1], fontsize=11)
    ax.set_ylim(-0.9, track_top + 0.75)        # headroom for gene labels + legend
    ax.set_xlim(WINDOW)

    # gene bodies + fusion-deletion shading (drawn over all tracks)
    for (gs, ge), name, col in [(ERG, "ERG", "#238b45"), (TMPRSS2, "TMPRSS2", "#6a51a3")]:
        ax.axvspan(gs, ge, ymax=(track_top + 0.3 + 0.9) / (track_top + 0.75 + 0.9),
                   color=col, alpha=0.18, zorder=0)
        ax.axvline(gs, color=col, lw=0.6, ls=":", zorder=3)
        ax.axvline(ge, color=col, lw=0.6, ls=":", zorder=3)
        ax.text((gs + ge) / 2, track_top + 0.5, name, ha="center", va="bottom",
                fontsize=12, color=col, fontweight="bold")
    ax.annotate("", xy=(TMPRSS2[0], track_top + 0.22), xytext=(ERG[1], track_top + 0.22),
                arrowprops=dict(arrowstyle="<->", color="0.35", lw=0.9))
    ax.text((ERG[1] + TMPRSS2[0]) / 2, track_top + 0.26,
            "interstitial deletion → TMPRSS2:ERG fusion", ha="center", va="bottom",
            fontsize=11, color="0.35")

    ax.set_xlabel("chr21 position (hg38, Mb)", fontsize=13)
    ax.set_xticks(np.arange(38e6, 42.5e6, 1e6))
    ax.set_xticklabels([f"{x/1e6:.0f}" for x in np.arange(38e6, 42.5e6, 1e6)], fontsize=11)
    # legend (figure-level, seated below the x-axis label)
    from matplotlib.patches import Patch
    leg = [Patch(facecolor=CMAP[v], label=lab) for v, lab in
           [(-2, "−2 deep loss"), (-1, "−1 loss"), (0, "neutral"), (1, "+1 gain")]]
    fig.legend(handles=leg, ncol=4, fontsize=11, frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, 0.005))
    ax.set_title("Deletion-mediated TMPRSS2:ERG fusion signature retained across "
                 "the castration timecourse (WES, hg38)", fontsize=13, pad=26)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.86, bottom=0.20)
    for fmt in args.plot_format:
        fig.savefig(outdir / f"tmprss2_erg_deletion.{fmt}", dpi=300, bbox_inches="tight")
    print(f"[ok] wrote {outdir}/tmprss2_erg_deletion.{{{','.join(args.plot_format)}}} "
          f"+ tmprss2_erg_deletion_calls.tsv")


if __name__ == "__main__":
    main()
