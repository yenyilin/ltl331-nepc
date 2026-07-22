#!/usr/bin/env python3
"""
compare_infercnv_vs_wes.py — side-by-side concordance of bulk WES CNV (from
nexus_cnv_loader.py) and scRNA-derived inferCNV calls (from infercnv_loader.py).

For each pair of matched samples (matched on the `timepoint` column in the
descriptors TSVs by default) the script:
  (a) bins both calls at --bin-size (default 1 Mbp; reuses the hg19/hg38 bin grid
      from plot_nexus_cnv_heatmap.py)
  (b) computes per-chromosome and genome-wide Pearson correlation
  (c) renders a stacked figure:
        - row 1: per-chrom Pearson r heatmap (rows = matched pairs, RdBu_r)
        - row 2: WES per-timepoint CNV heatmap
        - row 3: inferCNV per-timepoint CNV heatmap
      All three share the chromosome x-axis so chr1 in WES aligns with chr1
      in inferCNV and chr1's correlation in the top heatmap.

Outputs in --outdir:
  cnv_concordance.{pdf,png}            stacked 3-panel figure
  cnv_concordance_correlations.tsv     per-pair, per-chrom + genome-wide r

Example:
  python compare_infercnv_vs_wes.py \\
      --wes-calls-tsv data/nexus_calls.tsv \\
      --wes-descriptors-tsv data/nexus_descriptors.tsv \\
      --infercnv-calls-tsv data/infercnv_calls.tsv \\
      --infercnv-descriptors-tsv data/infercnv_descriptors.tsv \\
      --bin-size 1000000 --outdir nexus_plots --plot-format pdf png
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _utils import savefig, setup_style
from plot_nexus_cnv_heatmap import (
    CHROM_ORDER,
    DISCRETE_COLORS,
    HG19_SIZES,
    LEVELS,
    build_bins,
    paint_sample,
)


# --------------------------------------------------------------------------- #
# The concordance is restricted to the four same-study, newly-generated LTL331
# WES timepoints (GUNS-merged normal reference, same platform). The 8wk/12wk
# samples in the full descriptor are 927LT-normalized aCGH — a different
# platform + reference — and must NOT be paired into this figure.
SAME_STUDY_TIMEPOINTS = ["preCx", "16wk", "20wk", "22wk"]


def auto_pair(wes_descr, ic_descr, allowed_timepoints=None):
    """Match samples sharing a `timepoint`. Returns DataFrame
    [timepoint, wes_sample, infercnv_sample, n_pairs].

    If `allowed_timepoints` is given, only those timepoints are paired."""
    if "timepoint" not in wes_descr.columns or "timepoint" not in ic_descr.columns:
        raise SystemExit('both descriptors TSVs must have a "timepoint" column')

    allowed = set(allowed_timepoints) if allowed_timepoints else None
    pairs = []
    ic_by_tp = ic_descr.set_index("timepoint")["sample"].to_dict()
    seen_tp = set()
    for _, row in wes_descr.iterrows():
        tp = row.get("timepoint")
        if not isinstance(tp, str) or tp in seen_tp:
            continue
        if allowed is not None and tp not in allowed:
            continue
        if tp in ic_by_tp:
            pairs.append(
                {
                    "timepoint": tp,
                    "wes_sample": row["sample"],
                    "infercnv_sample": ic_by_tp[tp],
                }
            )
            seen_tp.add(tp)
    return pd.DataFrame(pairs)


def paint_to_matrix(calls, samples, bin_size, bins_df, chrom_offsets, chrom_nbins):
    """Returns (len(samples), n_bins) int8 matrix indexed by `samples`."""
    n = len(bins_df)
    mat = np.zeros((len(samples), n), dtype=np.int8)
    for i, s in enumerate(samples):
        sub = calls[calls["sample"] == s]
        if sub.empty:
            continue
        mat[i] = paint_sample(sub, bin_size, chrom_offsets, chrom_nbins, n)
    return mat


def per_chrom_corr(wes_row, ic_row, bins_df):
    """Per-chromosome + genome-wide Pearson r between two bin vectors.

    Pearson is undefined when one of the vectors has zero variance within the
    chromosome (e.g. both samples flat-neutral). Returns NaN there.
    """
    out = {}
    for chrom in CHROM_ORDER:
        cmask = (bins_df["chrom"] == chrom).values
        if not cmask.any():
            continue
        w, i = wes_row[cmask].astype(float), ic_row[cmask].astype(float)
        if w.std() == 0 or i.std() == 0:
            out[chrom] = float("nan")
        else:
            r, _ = stats.pearsonr(w, i)
            out[chrom] = r
    # global r across the full genome (more robust than averaging per-chrom)
    w, i = wes_row.astype(float), ic_row.astype(float)
    if w.std() == 0 or i.std() == 0:
        out["__genome__"] = float("nan")
    else:
        r, _ = stats.pearsonr(w, i)
        out["__genome__"] = r
    return out


# --------------------------------------------------------------------------- #
def _chrom_axis_decor(ax, bins_df):
    """Chromosome ticks on top + thin grey boundary lines."""
    boundaries, centers, labels = [], [], []
    for chrom, idx_series in bins_df.groupby("chrom", sort=False)["global_idx"]:
        idxs = idx_series.values
        boundaries.append(idxs[0] - 0.5)
        centers.append(idxs.mean())
        labels.append(chrom.replace("chr", ""))
    boundaries.append(bins_df["global_idx"].max() + 0.5)
    for b in boundaries[1:-1]:
        ax.axvline(b, color="#bbbbbb", lw=0.3)
    ax.set_xlim(boundaries[0], boundaries[-1])
    return boundaries, centers, labels


def _draw_cnv_heatmap(ax, matrix, row_labels, bins_df, panel_title):
    from matplotlib.colors import BoundaryNorm, ListedColormap

    cmap = ListedColormap([DISCRETE_COLORS[v] for v in LEVELS])
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=11)
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(panel_title, fontsize=9, loc="left", pad=4)
    _chrom_axis_decor(ax, bins_df)
    return im


def _draw_corr_heatmap(ax, corr_matrix, row_labels, bins_df):
    """corr_matrix shape: (n_pairs, n_bins) — but actually (n_pairs, n_chrom).
    We paint chromosome columns by repeating each chrom's r across its bins
    so the x-axis aligns with the CNV heatmaps below.
    """
    # corr_matrix is (n_pairs, n_chrom); broadcast to (n_pairs, n_bins)
    n_pairs, n_chrom = corr_matrix.shape
    n_bins = len(bins_df)
    bin_corr = np.full((n_pairs, n_bins), np.nan, dtype=float)
    chrom_to_col = {c: i for i, c in enumerate(CHROM_ORDER)}
    for chrom, idx_series in bins_df.groupby("chrom", sort=False)["global_idx"]:
        if chrom not in chrom_to_col:
            continue
        col = chrom_to_col[chrom]
        idxs = idx_series.values
        bin_corr[:, idxs] = corr_matrix[:, col][:, None]
    im = ax.imshow(
        bin_corr, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest"
    )
    ax.set_yticks(range(n_pairs))
    ax.set_yticklabels(row_labels, fontsize=11)
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(
        "per-chromosome Pearson r (WES vs inferCNV)", fontsize=9, loc="left", pad=4
    )
    # chromosome labels on top of the very top panel
    boundaries, centers, labels = _chrom_axis_decor(ax, bins_df)
    ax_top = ax.secondary_xaxis("top")
    ax_top.set_xticks(centers)
    ax_top.set_xticklabels(labels, fontsize=10)
    ax_top.tick_params(length=0)
    # annotate r values on the heatmap cells
    for i in range(n_pairs):
        for j, chrom in enumerate(CHROM_ORDER):
            r = corr_matrix[i, j]
            if not np.isfinite(r):
                continue
            x = centers[j]
            ax.text(
                x,
                i,
                f"{r:.2f}",
                ha="center",
                va="center",
                fontsize=6,
                color="white" if abs(r) > 0.55 else "#222",
            )
    return im


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__
    )
    p.add_argument("--wes-calls-tsv", required=True)
    p.add_argument("--wes-descriptors-tsv", required=True)
    p.add_argument("--infercnv-calls-tsv", required=True)
    p.add_argument("--infercnv-descriptors-tsv", required=True)
    p.add_argument("--bin-size", type=int, default=1_000_000)
    p.add_argument("--outdir", default="nexus_plots")
    p.add_argument(
        "--plot-format", nargs="+", default=["pdf"], choices=["pdf", "png", "svg"]
    )
    p.add_argument("--title", default=None)
    p.add_argument("--keep-sex-chrom", action="store_true",
                   help="keep chrX/chrY (default: drop them — inferCNV excludes the sex "
                        "chromosomes, so X/Y are always blank in the concordance)")
    p.add_argument("--timepoints", nargs="+", default=SAME_STUDY_TIMEPOINTS,
                   help='timepoints to pair (default: the four same-study WES timepoints '
                        f'{SAME_STUDY_TIMEPOINTS}). Pass "all" to pair every overlapping '
                        "timepoint (includes the 927LT-normalized 8wk/12wk aCGH samples).")
    args = p.parse_args()

    allowed_tp = None if args.timepoints == ["all"] else args.timepoints

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # inferCNV excludes the sex chromosomes → X/Y are always blank; drop them by default
    global CHROM_ORDER
    sizes = dict(HG19_SIZES)
    if not args.keep_sex_chrom:
        CHROM_ORDER = [c for c in CHROM_ORDER if c not in ("chrX", "chrY")]
        sizes = {k: v for k, v in HG19_SIZES.items() if k not in ("chrX", "chrY")}

    wes_calls = pd.read_csv(args.wes_calls_tsv, sep="\t")
    wes_descr = pd.read_csv(args.wes_descriptors_tsv, sep="\t")
    ic_calls = pd.read_csv(args.infercnv_calls_tsv, sep="\t")
    ic_descr = pd.read_csv(args.infercnv_descriptors_tsv, sep="\t")

    pairs = auto_pair(wes_descr, ic_descr, allowed_timepoints=allowed_tp)
    if pairs.empty:
        raise SystemExit("no overlapping timepoints between WES and inferCNV")
    # sort pairs by timepoint
    from _utils import timepoint_key

    pairs["__k"] = pairs["timepoint"].apply(timepoint_key)
    pairs = pairs.sort_values("__k").drop(columns="__k").reset_index(drop=True)
    print("Matched pairs:")
    print(pairs.to_string(index=False))

    bins_df = build_bins(args.bin_size, sizes)
    n_bins = len(bins_df)
    chrom_offsets, chrom_nbins = {}, {}
    for chrom, group in bins_df.groupby("chrom", sort=False):
        chrom_offsets[chrom] = int(group["global_idx"].min())
        chrom_nbins[chrom] = len(group)

    wes_mat = paint_to_matrix(
        wes_calls,
        pairs["wes_sample"].tolist(),
        args.bin_size,
        bins_df,
        chrom_offsets,
        chrom_nbins,
    )
    ic_mat = paint_to_matrix(
        ic_calls,
        pairs["infercnv_sample"].tolist(),
        args.bin_size,
        bins_df,
        chrom_offsets,
        chrom_nbins,
    )

    # per-chromosome + genome-wide correlation per pair
    corr_rows, corr_matrix = (
        [],
        np.full((len(pairs), len(CHROM_ORDER)), np.nan, dtype=float),
    )
    for i, (_, row) in enumerate(pairs.iterrows()):
        cd = per_chrom_corr(wes_mat[i], ic_mat[i], bins_df)
        for j, chrom in enumerate(CHROM_ORDER):
            corr_matrix[i, j] = cd.get(chrom, np.nan)
        rec = {
            "timepoint": row["timepoint"],
            "wes_sample": row["wes_sample"],
            "infercnv_sample": row["infercnv_sample"],
            "genome_pearson_r": cd["__genome__"],
        }
        rec.update({f"r_{c}": cd[c] for c in CHROM_ORDER if c in cd})
        corr_rows.append(rec)
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(outdir / "cnv_concordance_correlations.tsv", sep="\t", index=False)
    print("\nGenome-wide Pearson r per matched pair:")
    print(corr_df[["timepoint", "genome_pearson_r"]].to_string(index=False))

    # ----- 3-panel figure --------------------------------------------------- #
    # fig_w capped at 14in (not 20in): this figure is meant for a landscape
    # supplementary page (\linewidth ~9-10in there), and a narrower source
    # figure means less scale-down at print time -> larger effective font size
    # for the same absolute fontsize= values below (all vector, so text scales
    # with \includegraphics width like everything else in the PDF).
    plt = setup_style(font_size=9)
    n_pairs = len(pairs)
    fig_w = max(10.0, min(14.0, n_bins * 0.004 + 3.0))
    h_per_row = 0.32
    fig_h = (n_pairs * h_per_row * 3) + 3.0  # 3 panels each with n_pairs rows
    fig = plt.figure(figsize=(fig_w, fig_h))
    # reserve a dedicated colorbar column (col 1) so all three heatmaps (col 0) keep
    # IDENTICAL width and stay aligned. Attaching the colorbar to only the two CNV axes
    # (fig.colorbar(ax=[...])) would shrink just those two and break alignment with the
    # correlation panel above (sharex syncs data-limits, not physical axes width).
    gs = fig.add_gridspec(
        3,
        2,
        hspace=0.28,
        wspace=0.015,
        width_ratios=[1.0, 0.02],
        height_ratios=[
            max(n_pairs * 0.45, 2.0),
            max(n_pairs * 0.45, 2.0),
            max(n_pairs * 0.45, 2.0),
        ],
    )
    ax_corr = fig.add_subplot(gs[0, 0])
    ax_wes = fig.add_subplot(gs[1, 0], sharex=ax_corr)
    ax_ic = fig.add_subplot(gs[2, 0], sharex=ax_corr)
    cax = fig.add_subplot(gs[1:, 1])  # colorbar spans the two CNV rows, outside the panels

    pair_labels = [
        f"{r.timepoint}  r={cd:.2f}"
        for r, cd in zip(pairs.itertuples(), corr_df["genome_pearson_r"].values)
    ]
    # CNV-panel row labels: clean "LTL331_<stage>" (drop the raw "<sample> vs <ref>" names)
    stage_labels = [f"LTL331_{r.timepoint}" for r in pairs.itertuples()]

    _draw_corr_heatmap(ax_corr, corr_matrix, pair_labels, bins_df)
    im_wes = _draw_cnv_heatmap(
        ax_wes,
        wes_mat,
        stage_labels,
        bins_df,
        "Bulk WES (Nexus)",
    )
    im_ic = _draw_cnv_heatmap(
        ax_ic,
        ic_mat,
        stage_labels,
        bins_df,
        "scRNA inferCNV (HMMi6, Pnorm 0.5)",
    )

    # chromosome ruler along the BOTTOM of the inferCNV track too (the correlation panel
    # already labels chromosomes on top; a bottom copy gives a genomic reference under both
    # CNV tracks). secondary_xaxis is independent of the shared x, so it stays on this panel.
    _, centers, labels = _chrom_axis_decor(ax_ic, bins_df)
    ax_ic_bot = ax_ic.secondary_xaxis("bottom")
    ax_ic_bot.set_xticks(centers)
    ax_ic_bot.set_xticklabels(labels, fontsize=9)
    ax_ic_bot.tick_params(length=0)

    # shared discrete colorbar for the two CNV panels, in the reserved colorbar column
    cbar = fig.colorbar(im_wes, cax=cax, ticks=LEVELS)
    cbar.ax.set_yticklabels(
        ["big loss", "loss", "neutral", "gain", "high gain"], fontsize=8
    )
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)

    if args.title:
        fig.suptitle(args.title, fontsize=11, y=0.995)

    savefig(fig, outdir, "cnv_concordance", args.plot_format)
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
