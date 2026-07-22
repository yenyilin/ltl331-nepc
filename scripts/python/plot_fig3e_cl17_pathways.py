#!/usr/bin/env python3
"""
plot_fig3e_cl17_pathways.py — Fig 3E: representative enriched pathways in cluster 17
across three categories (EMT / Progenitor / Stemness), as a grouped horizontal NES bar
panel with FDR significance stars.

Consumes the C2:CGP GSEA output of per_cluster_gsea.py (run with --gmt c2.cgp):
  per_cluster_NES.tsv  (rows = pathway Term, cols = cluster)
  per_cluster_FDR.tsv  (same shape)
and a category->pathway mapping (data/fig3e_cl17_pathways.tsv; substring match). Plots the
NES of the focus cluster (default 17) for the matched pathways, grouped/coloured by category.

This is the cl17-specific companion to fig3d_nes_heatmap.py (the full 3C heatmap) — the
heatmap shows all clusters; this isolates cluster 17's representative programs for the caption.

Usage
  # 1) generate the C2:CGP NES first:
  python per_cluster_gsea.py --h5ad ltl331_base.h5ad --gmt c2.cgp.v7.5.1.symbols.gmt \\
      --use-raw --rank-mode vs_reference --reference-clusters 11 15 --outdir data/cgp_gsea
  # 2) plot 3E:
  python plot_fig3e_cl17_pathways.py --nes data/cgp_gsea/vs_reference/per_cluster_NES.tsv \\
      --fdr data/cgp_gsea/vs_reference/per_cluster_FDR.tsv --cluster 17 \\
      --pathways data/fig3e_cl17_pathways.tsv --outdir figures --stem fig3E
"""
import argparse
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402

CAT_COLOR = {"EMT": "#d73027", "Progenitor": "#fdae61", "Stemness": "#762a83"}


def load_map(path):
    cats = OrderedDict()
    for line in open(path):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        cats.setdefault(parts[0].strip(), []).append(parts[1].strip())
    return cats


def stars(fdr):
    if not np.isfinite(fdr):
        return ""
    return "***" if fdr < 1e-3 else "**" if fdr < 1e-2 else "*" if fdr < 0.05 else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nes", required=True, help="per_cluster_NES.tsv (Term x cluster)")
    ap.add_argument("--fdr", default=None, help="per_cluster_FDR.tsv (same shape)")
    ap.add_argument("--cluster", default="17")
    ap.add_argument("--pathways", default="data/fig3e_cl17_pathways.tsv")
    ap.add_argument("--max-per-cat", type=int, default=4, help="top |NES| per category")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--stem", default="fig3E")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    nes = pd.read_csv(args.nes, sep="\t", index_col=0)
    fdr = pd.read_csv(args.fdr, sep="\t", index_col=0) if args.fdr and os.path.exists(args.fdr) else None
    col = str(args.cluster)
    if col not in nes.columns:
        raise SystemExit(f"cluster {col!r} not a column of {args.nes}: {list(nes.columns)}")

    terms_lower = {t.lower(): t for t in nes.index.astype(str)}
    cats = load_map(args.pathways)
    rows = []
    for cat, patterns in cats.items():
        hits = []
        for pat in patterns:
            p = pat.lower()
            for tl, t in terms_lower.items():
                if p in tl:
                    v = float(nes.loc[t, col])
                    f = float(fdr.loc[t, col]) if fdr is not None and t in fdr.index else np.nan
                    hits.append((t, v, f))
        # dedup + top |NES|
        seen, uniq = set(), []
        for t, v, f in sorted(hits, key=lambda r: -abs(r[1])):
            if t not in seen:
                seen.add(t); uniq.append((t, v, f))
        if not uniq:
            print(f"[warn] category {cat!r}: no matching pathways in NES table")
        for t, v, f in uniq[:args.max_per_cat]:
            rows.append({"category": cat, "term": t, "NES": v, "FDR": f})
    if not rows:
        raise SystemExit("no pathways matched — check term names in --pathways vs the NES table")
    df = pd.DataFrame(rows)
    # order: category blocks, descending NES within each
    cat_order = [c for c in cats if c in set(df["category"])]
    df["cat_rank"] = df["category"].map({c: i for i, c in enumerate(cat_order)})
    df = df.sort_values(["cat_rank", "NES"], ascending=[True, True]).reset_index(drop=True)
    df.to_csv(os.path.join(args.outdir, f"{args.stem}_data.tsv"), sep="\t", index=False) \
        if os.path.isdir(args.outdir) else None

    plt = setup_style(font_size=7)
    os.makedirs(args.outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.0, 0.32 * len(df) + 1.6))
    y = np.arange(len(df))
    colors = [CAT_COLOR.get(c, "#666") for c in df["category"]]
    ax.barh(y, df["NES"], color=colors, edgecolor="black", lw=0.3, zorder=3)
    ax.axvline(0, color="0.5", lw=0.7)
    for yi, (_, r) in zip(y, df.iterrows()):
        s = stars(r["FDR"])
        if s:
            xoff = r["NES"] + (0.05 if r["NES"] >= 0 else -0.05)
            ax.text(xoff, yi, s, va="center", ha="left" if r["NES"] >= 0 else "right", fontsize=7)
    short = [t if len(t) <= 42 else t[:39] + "…" for t in df["term"]]
    ax.set_yticks(y); ax.set_yticklabels(short, fontsize=6)
    ax.set_xlabel(f"GSEA NES (cluster {col})")
    # category legend
    for c in cat_order:
        ax.barh([], [], color=CAT_COLOR.get(c, "#666"), label=c)
    ax.legend(loc="lower right", fontsize=7, frameon=False, title="category", title_fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Cluster {col}: representative enriched pathways", fontsize=8)

    savefig(fig, args.outdir, args.stem, args.plot_format)
    plt.close(fig)
    print(f"[ok] wrote {args.stem}.{{{','.join(args.plot_format)}}} + {args.stem}_data.tsv "
          f"to {args.outdir}/ ({len(df)} pathways across {len(cat_order)} categories)")


if __name__ == "__main__":
    main()
