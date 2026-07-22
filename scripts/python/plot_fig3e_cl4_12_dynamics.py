#!/usr/bin/env python3
"""
plot_fig3e_cl4_12_dynamics.py — Fig 3E: temporal dynamics of clusters 4 and 12.

Clusters 4 and 12 share a striking temporal arc: both expand pre-castration,
persist through the post-castration regression phase, then vanish as the NE
compartment takes over. The natural question is what role they play, survival
states that ride out castration, or precursors that seed the NE lineage. This
panel characterizes their dynamics to tell the two apart.

It shows each cluster's WITHIN-TIMEPOINT fraction across the 8 chronological
timepoints (preCx -> wk22R). Normalizing within timepoint is essential: timepoint
cell totals are very uneven (wk8 ~1.1k vs wk22R ~18.6k), so raw counts mislead.

What it shows: cluster 4 is a castration-SURVIVAL state, it peaks at the wk16 nadir
then is lost, NOT the NE origin (that is the cluster-17 bridge, shown separately in
Fig 4). The optional faint NE-aggregate reference shows the intermediates collapsing
AS the NE terminals take over, i.e. they "disappear into NE" but do not become it.

Source: per-cluster x timepoint counts (default: the DoubletFinder valid/post-QC table).

Example:
  python scripts/plot_fig3e_cl4_12_dynamics.py \
      --counts data/doubletfinder/per_cluster_timepoint_doublet_rates.tsv \
      --outdir figures/fig3
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

TP_ORDER = ["preCx", "wk2", "wk4", "wk8", "wk12", "wk16", "wk20R", "wk22R"]
TP_PRETTY = {
    "preCx": "preCx",
    "wk2": "2",
    "wk4": "4",
    "wk8": "8",
    "wk12": "12",
    "wk16": "16",
    "wk20R": "20R",
    "wk22R": "22R",
}
# NE terminal/lineage clusters (ASCL1- cl7/13/1/3/9/14 + ASCL1+ cl10) for the reference line
NE_CLUSTERS = [1, 3, 7, 9, 10, 13, 14]
CASTRATION_WINDOW = ("wk12", "wk16")  # weeks-12-16 castration window (regression nadir)
TERMINAL = ("wk20R", "wk22R")  # relapsed NE phase


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--counts", default="data/doubletfinder/per_cluster_timepoint_doublet_rates.tsv"
    )
    p.add_argument("--count-col", default="valid", help="count column (valid=post-QC)")
    p.add_argument("--targets", nargs="+", type=int, default=[4, 12])
    p.add_argument(
        "--show-ne",
        action="store_true",
        default=True,
        help='overlay faint NE-aggregate fraction as the "disappears in NE" reference',
    )
    p.add_argument("--no-ne", dest="show_ne", action="store_false")
    p.add_argument("--outdir", default="figures/fig3")
    p.add_argument("--stem", default="fig3E_cl4_12_dynamics")
    p.add_argument("--figwidth-mm", type=float, default=None,
                   help="author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa")
    p.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    args = p.parse_args()

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

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.counts, sep="\t")
    piv = (
        df.pivot_table(
            index="cluster", columns="timepoint", values=args.count_col, aggfunc="sum"
        )
        .reindex(columns=TP_ORDER)
        .fillna(0)
    )
    tot = piv.sum(0)
    frac = piv.div(tot, axis=1) * 100.0  # % of each timepoint
    ne_frac = frac.reindex([c for c in NE_CLUSTERS if c in frac.index]).sum(0)

    out = frac.loc[[c for c in args.targets if c in frac.index]].copy()
    out.loc["NE_aggregate"] = ne_frac
    out.to_csv(outdir / f"{args.stem}.tsv", sep="\t")
    print("=== within-timepoint % ===")
    print(out.round(1).to_string())

    from _utils import savefig, setup_style

    plt = setup_style()
    x = np.arange(len(TP_ORDER))
    colors = {4: "#d95f02", 12: "#7570b3"}
    titles = {4: "Cluster 4 (survival)", 12: "Cluster 12 (transient)"}
    n = len(args.targets)
    # Stack facets VERTICALLY (n rows x 1 col): each facet gets the full panel width so the
    # dual y-axes don't crush the plot, and the panel becomes ~square to pair with square 3C.
    fig, axes = plt.subplots(
        n, 1,
        figsize=figsize(3.4, 1.7 * n),
        squeeze=False, sharex=True,
    )
    for i, (ax, cl) in enumerate(zip(axes[:, 0], args.targets)):
        is_bottom = i == n - 1
        y = (
            frac.loc[cl, TP_ORDER].values
            if cl in frac.index
            else np.zeros(len(TP_ORDER))
        )
        bi = (TP_ORDER.index(CASTRATION_WINDOW[0]), TP_ORDER.index(CASTRATION_WINDOW[1]))
        ti = (TP_ORDER.index(TERMINAL[0]), TP_ORDER.index(TERMINAL[1]))
        ax.axvspan(bi[0] - 0.3, bi[1] + 0.3, color="#fff3cc", zorder=0)
        ax.axvspan(ti[0] - 0.3, ti[1] + 0.3, color="#e8f0e8", zorder=0)
        # target cluster on its OWN left axis (so cl12's small dynamics stay visible)
        ax.fill_between(x, 0, y, color=colors.get(cl, "#888"), alpha=0.25, zorder=2)
        (lc,) = ax.plot(
            x,
            y,
            "-o",
            color=colors.get(cl, "#888"),
            lw=1.6,
            ms=4,
            zorder=3,
            label=f"cluster {cl}",
        )
        ax.set_xticks(x)
        if is_bottom:  # sharex -> only the bottom facet carries the timepoint labels
            ax.set_xticklabels(
                [TP_PRETTY[t] for t in TP_ORDER], fontsize=6.5, rotation=45, ha="right",
                rotation_mode="anchor",
            )
            ax.set_xlabel("weeks post-castration", fontsize=9)
        ax.set_ylabel(
            f"cluster {cl} (% of tp)", fontsize=7, color=colors.get(cl, "#888")
        )
        ax.tick_params(axis="y", colors=colors.get(cl, "#888"), labelsize=6.5)
        ax.set_title(titles.get(cl, f"Cluster {cl}"), fontsize=8)
        ax.set_ylim(0, y.max() * 1.18 + 0.5)
        ax.spines[["top"]].set_visible(False)
        handles = [lc]
        # NE aggregate on a secondary right axis (0-100%) -> the "disappears in NE" reference
        if args.show_ne:
            ax2 = ax.twinx()
            (nl,) = ax2.plot(
                x,
                ne_frac[TP_ORDER].values,
                "--",
                color="0.55",
                lw=1.1,
                zorder=1,
                label="NE terminals (agg.)",
            )
            ax2.set_ylim(0, 105)
            ax2.set_ylabel("NE terminals (%)", fontsize=7, color="0.45")
            ax2.tick_params(axis="y", colors="0.45", labelsize=6.5)
            ax2.spines[["top"]].set_visible(False)
            handles.append(nl)
        if cl == args.targets[0]:
            top = ax.get_ylim()[1]
            ax.text(
                np.mean(bi),
                top * 0.97,
                "castration\nwindow",
                fontsize=6,
                color="#b8860b",
                ha="center",
                va="top",
            )
            ax.text(
                np.mean(ti),
                top * 0.97,
                "terminal\nNEPC",
                fontsize=6,
                color="#4a7a4a",
                ha="center",
                va="top",
            )
        ax.legend(handles=handles, frameon=False, fontsize=6.2, loc="upper left")
    fig.suptitle(
        "AR-low survival clusters expand during castration, then are lost in NEPC",
        fontsize=6.5,
        y = 0.9,#y=1.0,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    savefig(fig, outdir, args.stem, ("pdf", "png"))
    plt.close(fig)
    print(f"wrote {outdir}/{args.stem}.pdf/.png")


if __name__ == "__main__":
    main()
