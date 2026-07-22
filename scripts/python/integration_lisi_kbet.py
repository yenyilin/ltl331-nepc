#!/usr/bin/env python3
"""
integration_lisi_kbet.py — integration quality for the LTL331 + Gao + Li joint object.

Computes iLISI (batch mixing) + cLISI (cell-type preservation), and optionally kBET,
on the NAIVE-merge embedding (X_pca) vs the Harmony-corrected embedding
(X_harmony_dataset). The point is to show Harmony MIXES the three datasets without
OVER-CORRECTING away biology — i.e. iLISI improves while cLISI does not degrade.

Why we build our own label
--------------------------
cLISI needs a biological label that means the SAME thing across all three datasets:
  - the object's UCell scores + `pretty` are LTL331-only (NaN for Gao/Li),
  - `original_cluster` uses per-dataset (non-comparable) numbering,
so neither is usable directly. We assign a harmonized coarse class by scoring canonical
markers uniformly on .raw lognorm and taking the argmax:
    PRAD/luminal | NE | Basal/intermediate | Other(stromal/immune)
CAVEATS — this coarse argmax is a QC proxy, NOT a cell-type call:
    - It exists only to give cLISI a cross-dataset-comparable label for the batch-mixing
      metrics below. It is NOT the Fig 6A cell-type assignment: a hard argmax over-calls
      "Basal" (KRT5/TP63) on EMT-intermediate cells (~16% of LTL331), which contradicts the
      EMT-not-basal result. Fig 6A therefore uses the continuous adenocarcinoma->NE
      gradient, not these discrete labels. DO NOT read the "Basal" fraction here as biology.
    - The `Other` fraction is a ROUGH check of whether Gao/Li were restricted to tumor
      epithelium before integration. Because `Other` includes collagen/DCN (COL1A1, DCN, LUM),
      this argmax OVER-COUNTS non-epithelial content (ambient/low stromal signal is present in
      genuine epithelium). The definitive, DCN/collagen-corrected contamination metric lives in
      the clinical epithelial QC (Supplementary Figs. S19/S20), ~1.9%, not the ~17% a naive
      collagen-inclusive argmax implies. Cite S19/S20 for epithelial restriction, not this.

LISI directionality (RAW scale, as returned by harmonypy.compute_lisi)
  iLISI in [1, n_batch=3]   : HIGHER = better batch mixing
  cLISI in [1, n_class]     : LOWER  = better cell-type preservation (less mixing of types)
Good integration: iLISI(harmony) >> iLISI(pca)  AND  cLISI(harmony) ~= cLISI(pca).

Usage
  python integration_lisi_kbet.py --adata harmony.h5ad --out integration_qc/
  python integration_lisi_kbet.py --adata harmony.h5ad --label-key coarse_celltype  # reuse a label
  python integration_lisi_kbet.py --adata harmony.h5ad --kbet                        # also try kBET (needs scib + R kBET)
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc

# coarse, cross-dataset-comparable lineage markers (prostate adeno / NE / basal / non-epithelial)
MARKERS = {
    "PRAD":  ["AR", "KLK3", "KLK2", "FOLH1", "NKX3-1", "FOXA1", "KRT8", "KRT18"],
    "NE":    ["SYP", "CHGA", "CHGB", "NCAM1", "ENO2", "INSM1", "ASCL1", "SYT1", "SCG2"],
    "Basal": ["KRT5", "KRT14", "KRT15", "TP63"],
    "Other": ["PTPRC", "CD3D", "CD68", "VWF", "PECAM1", "COL1A1", "DCN", "LUM"],  # immune/endo/fibro
}


def assign_coarse(adata, use_raw=True):
    """Score each lineage on .raw lognorm and label each cell by the argmax program."""
    src = set(adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names)
    cols = []
    print("[label] scoring coarse lineage programs on", ".raw" if use_raw else ".X")
    for cls, genes in MARKERS.items():
        present = [g for g in genes if g in src]
        if not present:
            raise SystemExit(f"[abort] no markers found for {cls}; check var names / .raw")
        sc.tl.score_genes(adata, present, score_name=f"_sc_{cls}", use_raw=use_raw)
        cols.append(f"_sc_{cls}")
        print(f"        {cls:6s}: {len(present)}/{len(genes)} present -> {present}")
    S = adata.obs[cols].to_numpy()
    lab = np.array(list(MARKERS))[S.argmax(axis=1)]
    adata.obs["coarse_celltype"] = pd.Categorical(lab, categories=list(MARKERS))
    adata.obs.drop(columns=cols, inplace=True)
    return "coarse_celltype"


def lisi_medians(adata, rep, batch_key, label_key):
    """Median per-cell iLISI (batch) and cLISI (label) on one embedding (raw scale)."""
    from harmonypy import compute_lisi
    X = np.asarray(adata.obsm[rep])
    md = adata.obs[[batch_key, label_key]].astype(str)
    res = compute_lisi(X, md, [batch_key, label_key])  # (n_cells, 2)
    return float(np.median(res[:, 0])), float(np.median(res[:, 1]))


def try_kbet(adata, rep, batch_key, label_key):
    """kBET via scib (needs the R kBET pkg through rpy2); returns NaN + reason if unavailable."""
    try:
        import scib
        val = scib.metrics.kBET(adata, batch_key=batch_key, label_key=label_key,
                                type_="embed", embed=rep)
        return float(val), "ok"
    except Exception as e:  # scib missing, R kBET missing, or API drift
        return float("nan"), f"skipped: {type(e).__name__}: {e}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--adata", required=True)
    p.add_argument("--out", default="integration_qc")
    p.add_argument("--batch-key", default="dataset")
    p.add_argument("--label-key", default=None,
                   help="existing obs column to use as the cLISI cell-type label; "
                        "default = derive coarse_celltype from markers")
    p.add_argument("--reps", default="X_pca,X_harmony_dataset",
                   help="comma list: naive embedding first, corrected second")
    p.add_argument("--use-raw", action="store_true", default=True)
    p.add_argument("--no-raw", dest="use_raw", action="store_false")
    p.add_argument("--kbet", action="store_true", help="also attempt kBET (scib + R kBET)")
    p.add_argument("--save-h5ad", default=None, metavar="PATH",
                   help="write the object back WITH the derived coarse_celltype label "
                        "(so plot_umap_panels.py --color-by coarse_celltype can use it for Fig 6A)")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[load] {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"       {adata.n_obs:,} cells x {adata.n_vars:,} vars")

    # ---- gates ------------------------------------------------------------
    if not adata.obs_names.is_unique:
        raise SystemExit("[abort] obs_names not unique — fix before LISI (joins/graph assume unique).")
    reps = [r.strip() for r in args.reps.split(",") if r.strip()]
    for r in reps:
        if r not in adata.obsm:
            raise SystemExit(f"[abort] obsm[{r!r}] missing; have {list(adata.obsm)}")
    if args.batch_key not in adata.obs:
        raise SystemExit(f"[abort] batch key {args.batch_key!r} not in obs")
    print("[batch]", adata.obs[args.batch_key].value_counts().to_dict())

    # ---- cLISI label ------------------------------------------------------
    if args.label_key:
        if args.label_key not in adata.obs:
            raise SystemExit(f"[abort] --label-key {args.label_key!r} not in obs")
        if adata.obs[args.label_key].isna().any():
            raise SystemExit(f"[abort] label {args.label_key!r} has NaN — not usable across all cells")
        label_key = args.label_key
        print(f"[label] using existing obs[{label_key!r}]")
    else:
        label_key = assign_coarse(adata, use_raw=args.use_raw)

    # label QC: composition by dataset (+ vs standalone clustering if present)
    comp = pd.crosstab(adata.obs[label_key], adata.obs[args.batch_key])
    comp.to_csv(f"{args.out}/coarse_label_by_dataset.tsv", sep="\t")
    print("\n[label] coarse class x dataset:\n", comp)
    if "Other" in comp.index:
        for ds in comp.columns:
            frac = comp.loc["Other", ds] / comp[ds].sum()
            flag = "  <-- WARN: TME not removed?" if (ds != "LTL331" and frac > 0.25) else ""
            print(f"        Other fraction in {ds}: {frac:.1%}{flag}")
    if "original_cluster" in adata.obs:
        pd.crosstab(adata.obs["original_cluster"], adata.obs[label_key]).to_csv(
            f"{args.out}/coarse_label_vs_original_cluster.tsv", sep="\t")

    # persist the coarse lineage label for downstream plotting (Fig 6A)
    if args.save_h5ad:
        adata.write_h5ad(args.save_h5ad)
        print(f"[save] wrote object with obs[{label_key!r}] -> {args.save_h5ad}  "
              f"(use: plot_umap_panels.py --which categorical --color-by {label_key})")

    # ---- metrics: naive vs corrected -------------------------------------
    print(f"\n[lisi] batch={args.batch_key!r}  label={label_key!r}  (raw scale)")
    rows = []
    for rep in reps:
        il, cl = lisi_medians(adata, rep, args.batch_key, label_key)
        row = {"embedding": rep, "iLISI_median": round(il, 4), "cLISI_median": round(cl, 4)}
        if args.kbet:
            kb, why = try_kbet(adata, rep, args.batch_key, label_key)
            row["kBET"] = round(kb, 4) if kb == kb else kb
            if kb != kb:
                print(f"       [kBET {rep}] {why}")
        rows.append(row)
        print(f"       {rep:22s} iLISI={il:.3f}  cLISI={cl:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(f"{args.out}/lisi_kbet_metrics.tsv", sep="\t", index=False)

    # ---- verdict ----------------------------------------------------------
    nbatch = adata.obs[args.batch_key].nunique()
    verdict = {"n_batches": int(nbatch), "label_key": label_key,
               "metrics": rows, "scale": "raw (iLISI higher=better in [1,nbatch]; "
                                          "cLISI lower=better in [1,n_class])"}
    if len(rows) >= 2:
        naive, corr = rows[0], rows[1]
        d_i = corr["iLISI_median"] - naive["iLISI_median"]
        d_c = corr["cLISI_median"] - naive["cLISI_median"]
        verdict["delta_iLISI_corrected_minus_naive"] = round(d_i, 4)
        verdict["delta_cLISI_corrected_minus_naive"] = round(d_c, 4)
        good_mix = d_i > 0.1
        kept_bio = d_c < 0.3  # cLISI rose only modestly => types not smeared
        verdict["interpretation"] = (
            ("mixing IMPROVED" if good_mix else "mixing NOT improved")
            + "; " + ("biology PRESERVED" if kept_bio else "biology possibly OVER-CORRECTED")
        )
        print(f"\n[verdict] dILISI={d_i:+.3f}  dCLISI={d_c:+.3f}  -> {verdict['interpretation']}")
    with open(f"{args.out}/verdict.json", "w") as fh:
        json.dump(verdict, fh, indent=2)
    print(f"\n[ok] wrote {args.out}/lisi_kbet_metrics.tsv, coarse_label_by_dataset.tsv, verdict.json")


if __name__ == "__main__":
    main()
