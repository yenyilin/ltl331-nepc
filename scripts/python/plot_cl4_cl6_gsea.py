#!/usr/bin/env python3
"""
plot_cl4_cl6_gsea.py — prerank GSEA of clusters 4 and 6 vs. their trajectory
neighbours (Supplementary Fig. S10).

Four panels:
  A  cluster 4 vs cluster 0  (precursor)
  B  cluster 6 vs cluster 0  (precursor)
  C  cluster 4 vs cluster 17 (EMT successor)
  D  cluster 6 vs cluster 17 (EMT successor)

Method: prerank GSEA (gseapy.prerank). For each comparison we rank the FULL
transcriptome by the Wilcoxon test statistic (cluster-vs-reference, .raw lognorm)
and run GSEA against a gene-set library. Bars show normalized enrichment score
(NES, sign preserved); we report FDR q-value and leading-edge genes. This matches
the manuscript's "GSEA" wording and the original Fig 3D (fgsea) method family, and
— unlike ORA of a top-N list — recovers coordinate-but-diffuse programmes (e.g.
TGF-β, JAK/STAT) that no single cutoff would capture.

Gene-set library (--gene-sets): a local .gmt (offline; prerank is always offline).
Collections are inferred from MSigDB term-name PREFIXES, not the file (a merged gmt
carries no provenance):  HALLMARK_* -> hallmark ; GOBP_/GOCC_/GOMF_/GO_* -> go.
--collection {hallmark,go,both} filters which sets enter the test/plot. Use
`hallmark` for the cluster-4 signalling story (TNFα, p53, hypoxia, TGF-β, JAK/STAT),
`go` for the cluster-6 cilium story (GOCC motile cilium, axoneme), `both` for all.

Usage
  python plot_cl4_cl6_gsea.py \
      --h5ad ltl331_base.h5ad \
      --gene-sets refs/gmt/msigdb.hallmark_go.symbols.gmt \
      --collection both --out figures/supp/cl4_cl6_gsea --plot-format pdf png
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _utils import setup_style
    setup_style()
except Exception:
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})

# panel : (letter, test cluster, reference cluster, role)  -> 2x2 grid
COMPARISONS = [
    ("A", 4, 0, "precursor"),
    ("B", 6, 0, "precursor"),
    ("C", 4, 17, "EMT successor"),
    ("D", 6, 17, "EMT successor"),
]


# ---- gene-set handling (collection inferred from term-name prefix) ---------
def parse_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def collection_of(name):
    n = name.upper()
    if n.startswith("HALLMARK_"):
        return "hallmark"
    if n.startswith(("GOBP_", "GOCC_", "GOMF_", "GO_")):
        return "go"
    return "other"


def select_collection(sets, collection):
    if collection == "both":
        return dict(sets)
    keep = {k: v for k, v in sets.items() if collection_of(k) == collection}
    return keep


# ---- ranking + GSEA --------------------------------------------------------
def ranked_genes(adata, cluster_key, test, ref):
    """Full-transcriptome ranking by Wilcoxon score, test vs ref (.raw lognorm)."""
    n_all = adata.raw.n_vars if adata.raw is not None else adata.n_vars
    key = f"rgg_{test}_vs_{ref}"
    sc.tl.rank_genes_groups(
        adata, cluster_key, groups=[str(test)], reference=str(ref),
        method="wilcoxon", use_raw=True, n_genes=n_all, key_added=key)
    df = sc.get.rank_genes_groups_df(adata, group=str(test), key=key)
    df = df.dropna(subset=["names", "scores"]).drop_duplicates("names")
    rnk = (df[["names", "scores"]]
           .sort_values("scores", ascending=False)
           .reset_index(drop=True))
    return rnk


def run_gsea(rnk, gene_sets, min_size, max_size, seed):
    import gseapy as gp
    pre = gp.prerank(rnk=rnk, gene_sets=gene_sets, min_size=min_size,
                     max_size=max_size, permutation_num=1000, seed=seed,
                     threads=4, outdir=None, no_plot=True, verbose=False)
    res = pre.res2d.copy()
    ren = {}
    for c in res.columns:
        cl = str(c).lower()
        if cl == "term":
            ren[c] = "Term"
        elif cl.startswith("nes"):
            ren[c] = "NES"
        elif "fdr" in cl:
            ren[c] = "FDR"
        elif "lead" in cl:
            ren[c] = "Lead_genes"
    res = res.rename(columns=ren)
    # Term can arrive as "gene_sets__NAME" when a named gmt is used; strip it
    res["Term"] = res["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
    res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
    res["FDR"] = pd.to_numeric(res["FDR"], errors="coerce")
    return res.dropna(subset=["NES", "FDR"])


def clean_term(t, maxlen=44):
    t = re.sub(r"\s*\(GO:\d+\)$", "", str(t))
    t = re.sub(r"^(HALLMARK|GOBP|GOCC|GOMF|GO|REACTOME|KEGG)_", "", t, flags=re.I)
    t = t.replace("_", " ").strip()
    return t if len(t) <= maxlen else t[:maxlen - 1] + "…"


def print_terms(letter, test, ref, role, res, top_terms, fdr_max):
    hdr = f"Panel {letter}: cluster {test} vs cluster {ref} ({role})"
    print("\n" + "=" * len(hdr) + f"\n{hdr}\n" + "=" * len(hdr))
    sig = res[res["FDR"] <= fdr_max].sort_values("NES", ascending=False)
    if len(sig) == 0:
        print(f"  (no gene sets at FDR <= {fdr_max})")
        return
    for i, (_, r) in enumerate(sig.head(top_terms).iterrows(), 1):
        line = f"{i:2d}. {clean_term(r['Term'], 70)}  |  NES={r['NES']:+.2f}  |  FDR={r['FDR']:.1e}"
        print("   " + line)
        lg = r.get("Lead_genes", "")
        if isinstance(lg, str) and lg:
            print(f"        leading edge: {lg.replace(';', ', ')}")


def draw_panel(ax, letter, test, ref, role, res, top_terms, fdr_max, nes_max):
    ax.set_title(f"{letter}.  cluster {test} vs cluster {ref}  ({role})",
                 fontsize=9, fontweight="bold", loc="left")
    ax.set_xlabel("NES", fontsize=8)
    sig = res[res["FDR"] <= fdr_max].sort_values("NES", ascending=False)
    if len(sig) == 0:
        ax.text(0.5, 0.5, f"no sets at FDR ≤ {fdr_max}", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="0.5")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        return
    top = sig.head(top_terms).iloc[::-1]
    y = np.arange(len(top))
    nes = top["NES"].values
    norm = Normalize(vmin=-nes_max, vmax=nes_max)
    cols = plt.cm.coolwarm(norm(nes))
    ax.barh(y, nes, color=cols, edgecolor="0.3", linewidth=0.4)
    ax.axvline(0, color="0.3", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([clean_term(t) for t in top["Term"]], fontsize=7.5)
    # FDR annotation at the bar tip
    for yi, r in zip(y, top.itertuples()):
        v = r.NES
        ax.text(v + (0.06 * nes_max if v >= 0 else -0.06 * nes_max), yi,
                f"q={r.FDR:.0e}", va="center",
                ha="left" if v >= 0 else "right", fontsize=6.5, color="0.35")
    ax.set_xlim(-nes_max * 1.15, nes_max * 1.15)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--gene-sets", required=True, help="local .gmt path")
    ap.add_argument("--collection", choices=["hallmark", "go", "both"],
                    default="both")
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--out", default="figures/supp/cl4_cl6_gsea")
    ap.add_argument("--top-terms", type=int, default=12)
    ap.add_argument("--fdr-max", type=float, default=0.25)
    ap.add_argument("--min-size", type=int, default=10)
    ap.add_argument("--max-size", type=int, default=1500,
                    help="keep large GO sets (e.g. GOCC_CILIUM ~635 genes)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    ap.add_argument("--figwidth-mm", type=float, default=None,
                    help="author at this FINAL width in mm; default 154 mm (= \\textwidth). "
                         "Include with width=\\textwidth for 1:1 (true point-size) text; "
                         "the old 12-in canvas got halved when scaled to \\textwidth.")
    ap.add_argument("--figheight-mm", type=float, default=None)
    args = ap.parse_args()

    print("[version] plot_cl4_cl6_gsea rev7 (prerank GSEA + collection switch)")
    print(f"[load] {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    if adata.raw is None:
        print("[warn] .raw is None; ranking uses .X", file=sys.stderr)
    adata.obs[args.cluster_key] = (
        adata.obs[args.cluster_key].astype(str).astype("category"))
    present = set(adata.obs[args.cluster_key].cat.categories)
    needed = {str(c) for _, t, r, _ in COMPARISONS for c in (t, r)}
    missing = needed - present
    if missing:
        sys.exit(f"[error] clusters not found in {args.cluster_key}: {sorted(missing)}\n"
                 f"        present: {sorted(present)}")

    all_sets = parse_gmt(args.gene_sets)
    sets = select_collection(all_sets, args.collection)
    ncoll = {"hallmark": 0, "go": 0, "other": 0}
    for k in sets:
        ncoll[collection_of(k)] += 1
    print(f"[gene-sets] {args.gene_sets}: {len(all_sets)} total -> "
          f"collection='{args.collection}' keeps {len(sets)} "
          f"(hallmark={ncoll['hallmark']}, go={ncoll['go']}, other={ncoll['other']})")
    if not sets:
        sys.exit(f"[error] no sets match collection '{args.collection}'. "
                 f"Check term-name prefixes in {args.gene_sets}.")

    # first pass: compute all results, track global |NES| for a shared x-scale
    results = {}
    nes_max = 1.0
    for letter, test, ref, role in COMPARISONS:
        print(f"[{letter}] cluster {test} vs cluster {ref} ({role}) — ranking + GSEA")
        rnk = ranked_genes(adata, args.cluster_key, test, ref)
        res = run_gsea(rnk, sets, args.min_size, args.max_size, args.seed)
        results[letter] = (test, ref, role, res)
        sig = res[res["FDR"] <= args.fdr_max]
        if len(sig):
            nes_max = max(nes_max, float(np.nanmax(np.abs(sig["NES"].values))))
        print_terms(letter, test, ref, role, res, args.top_terms, args.fdr_max)

    def figsize(dw, dh):
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh:
            return (fw / 25.4, fh / 25.4)
        if fw:   # width only -> preserve authored aspect
            return (fw / 25.4, dh * (fw / 25.4) / dw)
        if fh:
            return (dw * (fh / 25.4) / dh, fh / 25.4)
        return (dw, dh)

    # authored at 6.06 in = 154 mm = \textwidth, so \includegraphics[width=\textwidth]
    # renders 1:1 (true point sizes) instead of shrinking the old 12-in canvas by ~0.5x.
    fig, axes = plt.subplots(2, 2, figsize=figsize(6.06, 5.6))
    axmap = {"A": axes[0, 0], "B": axes[0, 1], "C": axes[1, 0], "D": axes[1, 1]}
    all_res = []
    for letter, (test, ref, role, res) in results.items():
        draw_panel(axmap[letter], letter, test, ref, role, res,
                   args.top_terms, args.fdr_max, nes_max)
        all_res.append(res.assign(panel=letter, test=test, ref=ref,
                                  collection=[collection_of(t) for t in res["Term"]]))

    lib = os.path.basename(args.gene_sets)
    fig.suptitle(f"Prerank GSEA of clusters 4 and 6  (NES; {args.collection}; {lib})",
                 fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out = f"{args.out}_{args.collection}"
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    for fmt in args.plot_format:
        path = f"{out}.{fmt}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"[ok] wrote {path}")
    plt.close(fig)

    tsv = f"{out}_gsea.tsv"
    pd.concat(all_res).to_csv(tsv, sep="\t", index=False)
    print(f"[ok] wrote {tsv}")


if __name__ == "__main__":
    main()
