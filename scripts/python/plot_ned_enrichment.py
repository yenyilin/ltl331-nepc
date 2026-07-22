#!/usr/bin/env python3
"""
plot_ned_enrichment.py — Fig 6 panel C.

Renders the NED→NE enrichment from `validation_ne_enrichment.tsv`
(produced by `annotate_validation_phenotypes.py --report`): per clinical cohort, the
fraction of tumor epithelium called neuroendocrine (NE) in NED vs non-NED patients,
annotated with the Fisher odds ratio + p. The point: NE cells concentrate
in clinically NED-differentiated patients (OR >> 1).

Input columns (from the --report writer):
  dataset, NE_NED, NE_nonNED, rest_NED, rest_nonNED, odds_ratio, fisher_p

Usage:
  python plot_ned_enrichment.py \
      --enrichment data/integration_qc/validation_ne_enrichment.tsv \
      --out figures --stem fig6C_ned_enrichment --plot-format pdf png
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402

NED_COL, NONNED_COL = "#e7298a", "#9e9ac8"   # NE-pink (NED) vs muted purple (non-NED)


def _frac(num, den):
    den = num + den
    return 100.0 * num / den if den else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--enrichment", default="data/integration_qc/validation_ne_enrichment.tsv")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--stem", default="fig6C_ned_enrichment")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    df = pd.read_csv(args.enrichment, sep="\t")
    need = {"dataset", "NE_NED", "NE_nonNED", "rest_NED", "rest_nonNED", "odds_ratio", "fisher_p"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"[abort] {args.enrichment} missing columns {missing}; have {list(df.columns)}")
    df["ne_pct_NED"] = [_frac(r.NE_NED, r.rest_NED) for r in df.itertuples()]
    df["ne_pct_nonNED"] = [_frac(r.NE_nonNED, r.rest_nonNED) for r in df.itertuples()]

    plt = setup_style(font_size=8)
    fig, ax = plt.subplots(figsize=(max(3.2, 1.5 * len(df) + 1.2), 3.6))
    x = np.arange(len(df))
    w = 0.38
    ax.bar(x - w / 2, df["ne_pct_NED"], w, color=NED_COL, edgecolor="k", lw=0.4, label="NED")
    ax.bar(x + w / 2, df["ne_pct_nonNED"], w, color=NONNED_COL, edgecolor="k", lw=0.4, label="non-NED")

    ymax = float(np.nanmax([df["ne_pct_NED"].max(), df["ne_pct_nonNED"].max(), 1.0]))
    for i, r in df.iterrows():
        orr = r["odds_ratio"]
        or_str = "OR=∞" if not np.isfinite(orr) else (f"OR={orr:.0f}" if orr >= 10 else f"OR={orr:.1f}")
        ax.text(i, ymax * 1.04, f"{or_str}\np={r['fisher_p']}", ha="center", va="bottom", fontsize=6.5)
        # numeric labels on the bars
        ax.text(i - w / 2, df.loc[i, "ne_pct_NED"], f"{df.loc[i,'ne_pct_NED']:.0f}",
                ha="center", va="bottom", fontsize=6)
        ax.text(i + w / 2, df.loc[i, "ne_pct_nonNED"], f"{df.loc[i,'ne_pct_nonNED']:.1f}",
                ha="center", va="bottom", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(df["dataset"])
    ax.set_ylabel("NE epithelium (% of cohort's epithelium)")
    ax.set_ylim(0, ymax * 1.28)
    ax.set_title("Neuroendocrine epithelium is enriched in NED patients", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)

    os.makedirs(args.out, exist_ok=True)
    df[["dataset", "ne_pct_NED", "ne_pct_nonNED", "odds_ratio", "fisher_p"]].to_csv(
        os.path.join(args.out, f"{args.stem}.tsv"), sep="\t", index=False)
    savefig(fig, args.out, args.stem, args.plot_format)
    print(f"[ok] wrote {args.out}/{args.stem}.{{{','.join(args.plot_format)}}} (+ .tsv)")


if __name__ == "__main__":
    main()
