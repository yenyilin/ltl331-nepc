#!/usr/bin/env python3
"""
doublet_robustness.py — show, quantitatively AND visually, that the LTL331 results are not
driven by doublets. Consumes the unified tables from doublet_unify.py.

The argument has four legs, each a stat + a panel:
  (1) Doublets are rare and NOT enriched in the biologically load-bearing clusters.
      - per-cluster doublet rate vs the overall rate; cluster 17 (the EMT/bridge state)
        highlighted. Fisher exact per cluster (doublets-in-cluster vs rest), BH-corrected,
        to name which clusters ARE enriched (expected: the cycling engines) and confirm 17
        is not. Panel A.
  (2) Cluster composition is preserved after removal.
      - original vs valid per-cluster cell counts (Pearson/Spearman r, ~1) + cosine similarity
        of the composition vectors + max per-cluster fraction shift. Panel B (scatter).
  (3) The relative cluster structure barely shifts.
      - per-cluster fraction before vs after (paired). Panel C.
  (4) The temporal trajectory is unchanged — specifically cluster 17's wk12-16 peak survives.
      - cluster-17 fraction across timepoints, before vs after; per-timepoint composition
        correlation. Panel D.  Panel E: cluster x timepoint doublet-rate heatmap.

Outputs (OUT dir):
  doublet_robustness_stats.tsv     all numbers cited in the response
  cluster_enrichment.tsv           per-cluster Fisher OR + p + BH-q (doublet enrichment)
  doublet_robustness.{pdf,png}     5-panel supplementary figure

Usage:
  python doublet_robustness.py --unified data/doubletfinder/unified \
      --out data/doubletfinder/unified --target-cluster 17 \
      --palette data/cluster_colors_18.json --plot-format pdf png
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import savefig, setup_style  # noqa: E402


def bh(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--unified", default="data/doubletfinder/unified")
    ap.add_argument("--out", default=None)
    ap.add_argument("--target-cluster", type=int, default=17)
    ap.add_argument("--palette", default="data/cluster_colors_18.json")
    ap.add_argument(
        "--cluster-order",
        default="data/cluster_order.json",
        help="phenotypic cluster order (adeno->NE) for the A/C/E axes, to match "
        "every other figure; falls back to numeric sort if the file is absent",
    )
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()
    out = args.out or args.unified
    os.makedirs(out, exist_ok=True)
    tgt = args.target_cluster

    bc = pd.read_csv(f"{args.unified}/doublet_by_cluster.tsv", sep="\t").set_index(
        "cluster"
    )
    long = pd.read_csv(f"{args.unified}/doublet_counts_long.tsv", sep="\t")
    bc = bc.sort_index()
    # order the per-cluster axes (A/C/E) phenotypically (adeno->NE) like the rest of Fig 2/3,
    # instead of numeric 0..17. Order-independent stats (sums, Fisher, shifts) are unaffected.
    if args.cluster_order and os.path.exists(args.cluster_order):
        co = [int(c) for c in json.load(open(args.cluster_order)).get("order", [])]
        ordered = [c for c in co if c in bc.index] + [
            c for c in bc.index if c not in co
        ]
        bc = bc.reindex(ordered)
    clusters = list(bc.index)

    tot_o, tot_v = int(bc["original_n"].sum()), int(bc["valid_n"].sum())
    overall = 100 * (tot_o - tot_v) / tot_o

    # ---- (1) per-cluster enrichment (Fisher exact, doublet-in-cluster vs rest) ----
    try:
        from scipy.stats import fisher_exact

        have_scipy = True
    except Exception:
        have_scipy = False
    rows = []
    tot_d = tot_o - tot_v
    for c in clusters:
        d_c = max(int(bc.loc[c, "doublets"]), 0)  # clip negative (reassignment) to 0
        keep_c = int(bc.loc[c, "valid_n"])
        d_rest = tot_d - d_c
        keep_rest = (tot_v) - keep_c
        odds, p = (
            fisher_exact([[d_c, keep_c], [d_rest, keep_rest]])
            if have_scipy
            else (np.nan, np.nan)
        )
        rows.append(
            {
                "cluster": c,
                "doublet_pct": bc.loc[c, "doublet_pct"],
                "odds_ratio": round(float(odds), 3) if odds == odds else np.nan,
                "fisher_p": p,
            }
        )
    enr = pd.DataFrame(rows).set_index("cluster")
    if have_scipy:
        enr["bh_q"] = bh(enr["fisher_p"].values)
        enr["enriched"] = (enr["odds_ratio"] > 1) & (enr["bh_q"] < 0.05)
    enr["fisher_p"] = enr["fisher_p"].map(lambda x: f"{x:.2e}" if x == x else "NA")
    enr.reset_index().to_csv(f"{out}/cluster_enrichment.tsv", sep="\t", index=False)

    # ---- (2) composition preservation ----
    o_frac = bc["original_n"] / tot_o
    v_frac = bc["valid_n"] / tot_v
    pear = np.corrcoef(bc["original_n"], bc["valid_n"])[0, 1]
    sper = pd.Series(bc["original_n"]).rank().corr(pd.Series(bc["valid_n"]).rank())
    cos = float(
        np.dot(o_frac, v_frac) / (np.linalg.norm(o_frac) * np.linalg.norm(v_frac))
    )
    max_shift = float((100 * (v_frac - o_frac)).abs().max())
    max_shift_c = int((v_frac - o_frac).abs().idxmax())

    # ---- (4) temporal preservation ----
    piv_o = long.pivot_table(
        index="timepoint", columns="cluster", values="original_n", aggfunc="sum"
    ).fillna(0)
    piv_v = long.pivot_table(
        index="timepoint", columns="cluster", values="valid_n", aggfunc="sum"
    ).fillna(0)
    tps = long.sort_values("tp_order")["timepoint"].drop_duplicates().tolist()
    piv_o, piv_v = piv_o.reindex(tps), piv_v.reindex(tps)
    of = piv_o.div(piv_o.sum(1), axis=0)
    vf = piv_v.div(piv_v.sum(1), axis=0)
    per_tp_corr = {t: float(np.corrcoef(of.loc[t], vf.loc[t])[0, 1]) for t in tps}
    tgt_traj = pd.DataFrame({"original_frac": of[tgt], "valid_frac": vf[tgt]})

    # ---- write stats ----
    stats = {
        "overall_doublet_pct": round(overall, 2),
        "grand_original": tot_o,
        "grand_valid": tot_v,
        f"cluster{tgt}_doublet_pct": float(bc.loc[tgt, "doublet_pct"]),
        f"cluster{tgt}_original_n": int(bc.loc[tgt, "original_n"]),
        f"cluster{tgt}_valid_n": int(bc.loc[tgt, "valid_n"]),
        f"cluster{tgt}_enriched": (
            bool(enr.loc[tgt, "enriched"]) if "enriched" in enr else "scipy_absent"
        ),
        f"cluster{tgt}_OR": (
            float(enr.loc[tgt, "odds_ratio"])
            if enr.loc[tgt, "odds_ratio"] == enr.loc[tgt, "odds_ratio"]
            else None
        ),
        "composition_pearson_r": round(float(pear), 4),
        "composition_spearman_r": round(float(sper), 4),
        "composition_cosine": round(cos, 5),
        "max_cluster_frac_shift_pp": round(max_shift, 3),
        "max_shift_cluster": max_shift_c,
        "min_per_timepoint_composition_r": round(min(per_tp_corr.values()), 4),
        "per_timepoint_composition_r": {k: round(v, 4) for k, v in per_tp_corr.items()},
        "enriched_clusters": (
            [int(c) for c in enr.index[enr["enriched"]]]
            if "enriched" in enr
            else "scipy_absent"
        ),
    }
    pd.Series(stats, dtype=object).to_frame("value").to_csv(
        f"{out}/doublet_robustness_stats.tsv", sep="\t"
    )
    print(json.dumps(stats, indent=2, default=str))

    # ---- figure ----
    plt = setup_style()
    pal = {}
    if os.path.exists(args.palette):
        pal = {int(k): v for k, v in json.load(open(args.palette)).items()}
    col = [pal.get(c, "#bbbbbb") for c in clusters]
    hl = "#000000"  # target-cluster highlight: black, so it reads on ANY bar — incl. cl17,
    #                  whose own palette colour is magenta (#e7298a) and would hide a magenta outline

    # 3 rows: row1 = A (full width); row2 = B, C; row3 = D, E
    SUB_TITLE_FS = 15   # per-panel title size (large; canvas expands via savefig bbox='tight')
    SUP_TITLE_FS = 20   # main title size
    fig = plt.figure(figsize=(10, 11))
    gs = fig.add_gridspec(3, 2, hspace=0.80, wspace=0.30)

    # A: per-cluster doublet RATE  (doublets / that cluster's own cells)
    axA = fig.add_subplot(gs[0, :])
    bars = axA.bar(
        [str(c) for c in clusters], bc["doublet_pct"], color=col, edgecolor="k", lw=0.3
    )
    if tgt in clusters:
        bars[clusters.index(tgt)].set_edgecolor(hl)
        bars[clusters.index(tgt)].set_linewidth(1.8)
    axA.axhline(overall, ls="--", c="k", lw=0.8, label=f"overall {overall:.1f}%")
    if "enriched" in enr:
        for i, c in enumerate(clusters):
            if enr.loc[c, "enriched"]:
                axA.text(
                    i, bc.loc[c, "doublet_pct"] + 0.3, "*", ha="center", fontsize=9
                )
    axA.set_xlabel("cluster (adenocarcinoma → neuroendocrine)")
    axA.set_ylabel("doublet rate (%)\n= doublets ÷ cluster cells")
    axA.set_title(
        f"(A)  Per-cluster doublet RATE — within-cluster fraction called doublet\n"
        f"(cl{tgt} outlined; * = enriched, BH<0.05)",
        fontsize=SUB_TITLE_FS, fontweight="bold",
    )
    axA.legend(frameon=False, fontsize=7)

    # B: original vs valid counts (log-log)
    axB = fig.add_subplot(gs[1, 0])
    axB.scatter(
        bc["original_n"], bc["valid_n"], c=col, edgecolor="k", lw=0.3, s=28, zorder=3
    )
    if tgt in clusters:
        axB.scatter(
            [bc.loc[tgt, "original_n"]],
            [bc.loc[tgt, "valid_n"]],
            facecolor="none",
            edgecolor=hl,
            s=80,
            lw=1.8,
            zorder=4,
        )
    lim = [
        max(1, bc[["original_n", "valid_n"]].values.min() * 0.7),
        bc[["original_n", "valid_n"]].values.max() * 1.3,
    ]
    axB.plot(lim, lim, ls=":", c="grey", lw=0.8)
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(lim)
    axB.set_ylim(lim)
    axB.set_xlabel("original n")
    axB.set_ylabel("valid n")
    axB.set_title(f"(B)  Composition preserved\nr={pear:.3f}, cos={cos:.3f}",
                  fontsize=SUB_TITLE_FS, fontweight="bold")

    # C: per-cluster COMPOSITION (share of all cells) before vs after
    axC = fig.add_subplot(gs[1, 1])
    x = np.arange(len(clusters))
    w = 0.4
    axC.bar(
        x - w / 2,
        100 * o_frac.values,
        w,
        label="original",
        color="#9ecae1",
        edgecolor="k",
        lw=0.3,
    )
    axC.bar(
        x + w / 2,
        100 * v_frac.values,
        w,
        label="valid",
        color="#fc9272",
        edgecolor="k",
        lw=0.3,
    )
    axC.set_xticks(x)
    axC.set_xticklabels([str(c) for c in clusters], fontsize=6)
    axC.set_xlabel("cluster (adenocarcinoma → neuroendocrine)")
    axC.set_ylabel("% of all cells\n= cluster cells ÷ dataset")
    axC.set_title(
        f"(C)  Cluster COMPOSITION (share of dataset) \nbefore vs after —"
        f"max shift {max_shift:.2f} pp @ cl{max_shift_c}",
        fontsize=SUB_TITLE_FS, fontweight="bold",
    )
    axC.legend(frameon=False, fontsize=7)

    # D: target-cluster temporal trajectory. Annotate the per-timepoint cell count and
    # ring the timepoints where the cluster has ZERO cells, so y=0 reads as "absent"
    # (e.g. cl17: preCx/wk2/wk4/wk22) rather than "tiny fraction".
    axD = fig.add_subplot(gs[2, 0])
    cnt_o = (piv_o[tgt].reindex(tps).fillna(0).astype(int)
             if tgt in piv_o.columns else pd.Series(0, index=tps)).to_numpy()
    yo = 100 * tgt_traj["original_frac"].to_numpy()
    axD.plot(tps, yo, "-o", ms=3, c="#3182bd", label="original")
    axD.plot(tps, 100 * tgt_traj["valid_frac"].to_numpy(), "--s", ms=3, c=hl, label="valid")
    absent = cnt_o == 0
    if absent.any():
        axD.scatter(np.asarray(tps, dtype=object)[absent], yo[absent], s=55,
                    facecolors="none", edgecolors="#d62728", linewidths=1.3, zorder=5,
                    label="0 cells (cluster absent)")
    for i, t in enumerate(tps):
        axD.annotate(f"n={cnt_o[i]}", (t, yo[i]), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=6,
                     color=("#d62728" if absent[i] else "0.35"),
                     fontweight=("bold" if absent[i] else "normal"))
    axD.set_xlabel("timepoint")
    axD.set_ylabel(f"cluster {tgt}  (% of timepoint)")
    axD.set_title(f"(D)  Cluster {tgt} trajectory preserved",
                  fontsize=SUB_TITLE_FS, fontweight="bold")
    axD.tick_params(axis="x", rotation=45)
    axD.margins(y=0.18)   # headroom for the n= labels
    axD.legend(frameon=False, fontsize=7)

    # E: cluster x timepoint doublet-rate heatmap (rows in the same phenotypic order as A/C)
    axE = fig.add_subplot(gs[2, 1])
    hm = long.pivot_table(
        index="cluster", columns="timepoint", values="doublet_pct", aggfunc="mean"
    ).reindex(columns=tps)
    hm = hm.reindex(index=[c for c in clusters if c in hm.index][::-1])
    im = axE.imshow(
        hm.values,
        aspect="auto",
        cmap="Reds",
        vmin=0,
        vmax=np.nanpercentile(hm.values, 95),
    )
    axE.set_yticks(range(len(hm.index)))
    axE.set_yticklabels(hm.index, fontsize=6)
    axE.set_xticks(range(len(tps)))
    axE.set_xticklabels(tps, rotation=45, fontsize=6)
    axE.set_title("(E)  doublet % (cluster x timepoint)",
                  fontsize=SUB_TITLE_FS, fontweight="bold")
    fig.colorbar(im, ax=axE, fraction=0.046, pad=0.04).set_label(
        "doublet %", fontsize=7
    )

    fig.suptitle(
        "Cluster structure and trajectory are robust to doublet exclusion:\n"
        "cluster 17 is not doublet-enriched",
        fontsize=SUP_TITLE_FS,
        fontweight="bold",
        y=1.04,
    )
    savefig(fig, out, "doublet_robustness", args.plot_format)
    print(
        f"\n[ok] wrote doublet_robustness_stats.tsv, cluster_enrichment.tsv, "
        f"doublet_robustness.{{{','.join(args.plot_format)}}} to {out}/"
    )


if __name__ == "__main__":
    main()
