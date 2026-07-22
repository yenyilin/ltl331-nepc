#!/usr/bin/env python3
"""
plot_cl17_stemness_programs.py — Supplementary "cluster 17 stemness" figure, PANEL A
(per-cluster program scores). Redesign of the muddled original Supp Fig 3 B/C.

WHAT IT SHOWS
  A compact heatmap: rows = five curated de-differentiation programs (EMT, prostate cancer
  stem cell, luminal progenitor, neural-crest stem cell, mesenchymal stem cell), columns =
  clusters in the canonical adeno->NE order. Color = per-cluster mean program score, z-scored
  across clusters per row so the cluster carrying each program stands out. Cluster 17 (the
  EMT/stem bridge) is boxed and carries significance stars from a cluster-17-vs-rest test.

WHY A SCORE PANEL (not the original per-gene dot plots for this question)
  The individual "stem" markers are non-specific (integrins, CD44, ...). The defensible claim
  is PROGRAM-LEVEL co-activation with a statistic, not single genes. The gene-level receipts
  live in the SEPARATE dot plot (data/cl17_stemness_markers.tsv -> plot_marker_dotplot.py),
  which keeps "program enriched" and "gene expressed" as two clean, non-conflicting messages.

SCORING: sc.tl.score_genes on .raw lognorm (LTL331 .X is scaled — see [[ltl331_h5ad_schema]]).
STATS:  per program, Mann-Whitney U of focus-cluster cells vs all other cells, Benjamini-
        Hochberg-adjusted across the five programs (q in the stats TSV; stars on the heatmap).

Outputs (in --out): <stem>.{pdf,png}, <stem>_per_cluster_scores.tsv (raw means + row z),
<stem>_cl17_stats.tsv (MWU U/p/q, focus vs rest medians, n genes found/total per program).

Usage
  python plot_cl17_stemness_programs.py --h5ad ltl331_base.h5ad \\
      --signatures data/cl17_stemness_signatures.tsv --focus-cluster 17 \\
      --cluster-order data/cluster_order.json --out figures/figS8_stemness --stem figS8A_programs

Panel B (the gene-level dot plot) is rendered by the canonical Fig-2C renderer:
  python plot_marker_dotplot.py --h5ad ltl331_base.h5ad --markers data/cl17_stemness_markers.tsv \\
      --cluster-order data/cluster_order.json --out figures/figS8_stemness --stem figS8B_markers
"""
import argparse
import json
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig, numkey, to_dense_1d  # noqa: E402


def load_sigs(path):
    """Read a <program>\\t<gene> TSV (preserving program order, # comments ok)."""
    sigs = OrderedDict()
    for line in open(path):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        sigs.setdefault(parts[0].strip(), []).append(parts[1].strip())
    return sigs


def bh_adjust(pvals):
    """Benjamini-Hochberg q-values (no statsmodels dep)."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def stars(q):
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 5e-2 else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--signatures", default="data/cl17_stemness_signatures.tsv")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--focus-cluster", default="17")
    ap.add_argument("--cluster-order", default=None)
    ap.add_argument("--exclude-clusters", nargs="*", default=["18"])
    ap.add_argument("--out", default="figures/figS8_stemness")
    ap.add_argument("--stem", default="figS8A_programs")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    plt = setup_style()

    adata = sc.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    if args.exclude_clusters:
        adata = adata[~adata.obs[args.cluster_key].isin(args.exclude_clusters)].copy()

    # cluster order: canonical file, intersected with what is present; else numeric
    present = list(pd.unique(adata.obs[args.cluster_key]))
    if args.cluster_order:
        co = [str(c) for c in json.load(open(args.cluster_order)).get("order", [])]
        clusters = [c for c in co if c in present] + sorted(
            [c for c in present if c not in co], key=numkey)
    else:
        clusters = sorted(present, key=numkey)

    # score programs on raw lognorm; report missing genes
    sigs = load_sigs(args.signatures)
    use_raw = adata.raw is not None
    src = set(map(str, (adata.raw.var_names if use_raw else adata.var_names)))
    found = OrderedDict()
    for name, genes in sigs.items():
        gpresent = [g for g in genes if g in src]
        missing = [g for g in genes if g not in src]
        if missing:
            print(f"  [{name}] missing {len(missing)}/{len(genes)}: {','.join(missing)}")
        if len(gpresent) < 2:
            print(f"  [{name}] SKIPPED (<2 genes found)")
            continue
        sc.tl.score_genes(adata, gpresent, score_name=name, use_raw=use_raw)
        found[name] = (len(gpresent), len(genes))
    programs = list(found.keys())

    cl = adata.obs[args.cluster_key].values
    focus = args.focus_cluster

    # per-cluster mean score matrix (programs x clusters)
    M = np.zeros((len(programs), len(clusters)))
    for i, pname in enumerate(programs):
        s = adata.obs[pname].values
        for j, c in enumerate(clusters):
            M[i, j] = float(np.mean(s[cl == c])) if np.any(cl == c) else np.nan

    # row z-score across clusters (for color contrast only)
    Z = (M - np.nanmean(M, axis=1, keepdims=True)) / (np.nanstd(M, axis=1, keepdims=True) + 1e-9)

    # cluster-17-vs-rest stats per program
    rows = []
    pvals = []
    for pname in programs:
        s = adata.obs[pname].values
        a, b = s[cl == focus], s[cl != focus]
        U, p = mannwhitneyu(a, b, alternative="greater")
        pvals.append(p)
        rows.append(dict(program=pname, genes_found=found[pname][0], genes_total=found[pname][1],
                         focus_median=float(np.median(a)), rest_median=float(np.median(b)),
                         U=float(U), p_greater=float(p)))
    qvals = bh_adjust(pvals)
    for r, q in zip(rows, qvals):
        r["q_BH"] = float(q)
    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(os.path.join(args.out, f"{args.stem}_cl17_stats.tsv"), sep="\t", index=False)

    # raw means + z to TSV
    sc_df = pd.DataFrame(M, index=programs, columns=clusters)
    z_df = pd.DataFrame(Z, index=[f"{p}__z" for p in programs], columns=clusters)
    pd.concat([sc_df, z_df]).to_csv(
        os.path.join(args.out, f"{args.stem}_per_cluster_scores.tsv"), sep="\t")

    # ---- heatmap ----
    fig_w = max(5.0, 0.34 * len(clusters) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, 0.42 * len(programs) + 1.6))
    vmax = np.nanmax(np.abs(Z))
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, fontsize=7)
    ax.set_yticks(range(len(programs)))
    ax.set_yticklabels([p.replace("_", " ") for p in programs], fontsize=8)
    ax.set_xlabel("cluster (adenocarcinoma → neuroendocrine axis)", fontsize=8)

    # box + stars on the focus column
    if focus in clusters:
        fj = clusters.index(focus)
        ax.add_patch(plt.Rectangle((fj - 0.5, -0.5), 1, len(programs),
                                   fill=False, ec="black", lw=1.6))
        for r, q in zip(rows, qvals):
            i = programs.index(r["program"])
            st = stars(q)
            if st:
                ax.text(fj, i, st, ha="center", va="center", fontsize=8, fontweight="bold",
                        color="black")
        ax.get_xticklabels()[fj].set_fontweight("bold")

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("row z-scored mean program score", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Cluster 17 co-activates EMT / stem / progenitor programs\n"
                 f"(stars = cl{focus} vs rest, Mann–Whitney, BH-adj)", fontsize=8)

    fig.tight_layout()
    savefig(fig, args.out, args.stem, args.plot_format)
    print(f"  programs scored: {len(programs)} | clusters: {len(clusters)} | focus cl{focus}")
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
