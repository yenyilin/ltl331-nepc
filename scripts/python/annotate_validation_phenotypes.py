#!/usr/bin/env python3
"""
annotate_validation_phenotypes.py — add patient-level phenotype covariates to the joint
LTL331 + Gao + Li object, keyed by `orig.ident`.

Adds obs columns:
  dataset            LTL331 | Gao | Li         (derived from orig.ident if not present)
  patient_phenotype  detailed per-patient label (e.g. Gao_NEPC_smallcell, Li_mCRPC_pureSmallCell)
  ned_status         NED | non-NED | PDX       (neuroendocrine-differentiation status)
  hormone_status     HSPC | CRPC | PDX_castration_series
  path_note          free-text pathology note (provenance for the call)

Ground truth verified 2026-06-16 against the source papers:

  Gao = Dong et al., Comms Biol 2020 — 6 CRPC needle biopsies. Per-biopsy pathology
        (Results / Fig 1B / Table 1; NE scoring Fig 1C):
          #1,#4 = classical PCA (PRAD);  #3 = PIN (paper flags as biopsy-sampling artifact);
          #2,#5,#6 = small-cell carcinoma morphology + NED (NEPC).  All 6 clinically CRPC.
  Li  = Wang et al., iScience 2022 — used ONLY the 4 native patients (GSA HRA002145),
        object IDs Li_1/2/3/6 (Li reused Gao P2/P5/P6 as their P7/P4/P5 → excluded, no
        double-counting; Li_3 Drop-seq):
          Li_1=P1 mHSPC adeno; Li_2=P2 nmHSPC mixed-NEPC; Li_3=P3 mCRPC mixed-NEPC;
          Li_6=P6 mCRPC pure small-cell NEPC.

NOTE on keys: assumes Gao object idents are Gao_<N> == Dong patient #N, and Li_<N> ==
Wang patient P<N>. The script PRINTS the orig.ident -> mapping table and WARNS on any
ident it cannot place, so a numbering surprise is caught, not silently mislabeled.

Usage:
  python annotate_validation_phenotypes.py --adata harmony.h5ad --dry-run   # report only
  python annotate_validation_phenotypes.py --adata harmony.h5ad --out harmony.annot.h5ad
"""
import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc

# patient-number -> phenotype, keyed within each dataset (see docstring for sources)
GAO = {
    1: dict(phenotype="Gao_PRAD",           ned="non-NED", hormone="CRPC", note="classical PCA (Fig 1B)"),
    2: dict(phenotype="Gao_NEPC_smallcell", ned="NED",     hormone="CRPC", note="small-cell carcinoma; majority-NE epithelium (Fig 1C)"),
    3: dict(phenotype="Gao_PIN",            ned="non-NED", hormone="CRPC", note="PIN — paper notes likely biopsy-sampling inaccuracy"),
    4: dict(phenotype="Gao_PRAD",           ned="non-NED", hormone="CRPC", note="classical PCA (Fig 1B)"),
    5: dict(phenotype="Gao_NEPC_smallcell", ned="NED",     hormone="CRPC", note="small-cell carcinoma; majority-NE epithelium (Fig 1C)"),
    6: dict(phenotype="Gao_NEPC_smallcell", ned="NED",     hormone="CRPC", note="small-cell carcinoma; partial-NE epithelium (Fig 1C)"),
}
LI = {
    1: dict(phenotype="Li_mHSPC_adeno",         ned="non-NED", hormone="HSPC", note="mHSPC adenocarcinoma (Table S1)"),
    2: dict(phenotype="Li_nmHSPC_mixedNEPC",    ned="NED",     hormone="HSPC", note="nmHSPC mixed-NEPC; Syn+/CD56+/NSE+ (Table S1)"),
    3: dict(phenotype="Li_mCRPC_mixedNEPC",     ned="NED",     hormone="CRPC", note="mCRPC mixed-NEPC; CgA+ (Table S1)"),
    6: dict(phenotype="Li_mCRPC_pureSmallCell", ned="NED",     hormone="CRPC", note="mCRPC pure small-cell NEPC (Table S1)"),
}


def ensure_coarse(adata, use_raw=True):
    """Return the coarse-lineage label column, computing it if absent.

    Reuses assign_coarse() from integration_lisi_kbet.py so the PRAD/NE/Basal/Other
    marker definitions are a single source of truth (no duplicated marker lists)."""
    if "coarse_celltype" in adata.obs:
        return "coarse_celltype"
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from integration_lisi_kbet import assign_coarse  # noqa: E402  (sibling script)
    if use_raw and adata.raw is None:
        print("[report] no .raw — scoring coarse lineage on .X (check it is lognorm!)")
    return assign_coarse(adata, use_raw=(use_raw and adata.raw is not None))


def report_validation(adata, dataset_key, out_dir, use_raw=True):
    """Emit coarse_celltype x ned_status PER clinical dataset — the
    validation claim (NED patients' epithelium scores NE; non-NED scores luminal).

    Prints column-normalized tables (each ned_status sums to 1) and writes one tidy
    long-form TSV (dataset, ned_status, coarse_celltype, n, frac_within_ned) for plotting."""
    label = ensure_coarse(adata, use_raw=use_raw)
    ned = "ned_status"
    rows, cts = [], {}
    print("\n[report] coarse lineage x NED status, per dataset "
          "(column % = composition within each NED group):")
    datasets = sorted(adata.obs[dataset_key].astype(str).unique())
    for ds in datasets:
        mask = adata.obs[dataset_key].astype(str) == str(ds)
        ct = pd.crosstab(adata.obs.loc[mask, label], adata.obs.loc[mask, ned])
        if ct.values.sum() == 0:
            continue
        cts[ds] = ct
        frac = ct.div(ct.sum(axis=0).replace(0, pd.NA), axis=1)
        print(f"\n  === {ds} (n={int(mask.sum()):,}) ===")
        print((frac * 100).round(1).fillna(0).to_string())
        for coarse in ct.index:
            for nd in ct.columns:
                n = int(ct.loc[coarse, nd])
                rows.append(dict(dataset=ds, ned_status=nd, coarse_celltype=coarse, n=n,
                                 frac_within_ned=(None if ct[nd].sum() == 0
                                                  else round(n / ct[nd].sum(), 4))))
    tidy = pd.DataFrame(rows)
    path = os.path.join(out_dir, "validation_coarse_by_ned.tsv")
    os.makedirs(out_dir, exist_ok=True)
    tidy.to_csv(path, sep="\t", index=False)
    print(f"\n[report] wrote {path}  ({len(tidy)} rows; plot-ready long form)")

    _ne_enrichment(cts, out_dir)
    return label


def _ne_enrichment(cts, out_dir):
    """NE-vs-rest x (NED vs non-NED) 2x2 odds ratio + Fisher exact p, per clinical cohort.
    Quantifies the claim 'NE epithelium tracks the NED pathology call'. Skips cohorts
    lacking both NED and non-NED (e.g. the LTL331 PDX, which is PDX-status only)."""
    try:
        from scipy.stats import fisher_exact
    except Exception as e:
        print(f"[report] scipy unavailable ({type(e).__name__}) — skipping NE-enrichment test")
        return
    rows = []
    for ds, ct in cts.items():
        if "NE" not in ct.index or not {"NED", "non-NED"}.issubset(set(ct.columns)):
            continue
        ne_ned, ne_non = int(ct.loc["NE", "NED"]), int(ct.loc["NE", "non-NED"])
        rest_ned = int(ct["NED"].sum()) - ne_ned
        rest_non = int(ct["non-NED"].sum()) - ne_non
        # table: [[NE&NED, NE&nonNED], [rest&NED, rest&nonNED]]
        odds, p = fisher_exact([[ne_ned, ne_non], [rest_ned, rest_non]])
        rows.append(dict(dataset=ds, NE_NED=ne_ned, NE_nonNED=ne_non,
                         rest_NED=rest_ned, rest_nonNED=rest_non,
                         odds_ratio=round(float(odds), 2), fisher_p=f"{p:.3e}"))
    if not rows:
        print("[report] no cohort has both NED and non-NED — no enrichment test")
        return
    df = pd.DataFrame(rows)
    path = os.path.join(out_dir, "validation_ne_enrichment.tsv")
    df.to_csv(path, sep="\t", index=False)
    print("\n[report] NE enrichment in NED epithelium (NE-vs-rest x NED-vs-nonNED):")
    print(df.to_string(index=False))
    print(f"[report] wrote {path}")


def classify_dataset(ident):
    s = str(ident)
    if s.lower().startswith("gao") or "patient" in s.lower():
        return "Gao"
    if s.lower().startswith("li"):
        return "Li"
    if re.fullmatch(r"\d+", s):       # bare timepoint number => LTL331
        return "LTL331"
    return "UNKNOWN"


def patient_number(ident):
    m = re.search(r"(\d+)\s*$", str(ident))   # trailing integer: Gao_3, "patient #3", Li_6
    return int(m.group(1)) if m else None


def resolve(ident):
    """Return (dataset, dict) for one orig.ident, or (dataset, None) if unmappable."""
    ds = classify_dataset(ident)
    if ds == "LTL331":
        return ds, dict(phenotype=f"LTL331_{ident}", ned="PDX",
                        hormone="PDX_castration_series", note="PDX trajectory timepoint")
    if ds == "Gao":
        return ds, GAO.get(patient_number(ident))
    if ds == "Li":
        return ds, LI.get(patient_number(ident))
    return ds, None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--adata", required=True)
    p.add_argument("--out", default=None, help="output h5ad (default: <in>.annot.h5ad)")
    p.add_argument("--orig-key", default="orig.ident")
    p.add_argument("--dataset-key", default="dataset")
    p.add_argument("--dry-run", action="store_true", help="report mapping, do not write")
    p.add_argument("--report", action="store_true",
                   help="also emit coarse_celltype x ned_status per dataset"
                        "; computes coarse_celltype from markers if absent")
    p.add_argument("--report-dir", default="data/integration_qc",
                   help="where to write validation_coarse_by_ned.tsv (with --report)")
    p.add_argument("--no-raw", dest="use_raw", action="store_false", default=True,
                   help="score coarse lineage on .X instead of .raw (with --report)")
    args = p.parse_args()

    print(f"[load] {args.adata}")
    adata = sc.read_h5ad(args.adata)
    if args.orig_key not in adata.obs:
        raise SystemExit(f"[abort] obs[{args.orig_key!r}] not found; have {list(adata.obs)[:10]}...")

    idents = adata.obs[args.orig_key].astype(str)
    uniq = sorted(idents.unique())

    # build the per-ident mapping table and surface it for inspection
    rows, unmapped = [], []
    mapping = {}
    for ident in uniq:
        ds, d = resolve(ident)
        n = int((idents == ident).sum())
        if d is None:
            unmapped.append(ident)
            rows.append(dict(orig_ident=ident, dataset=ds, n=n, patient_phenotype="<<UNMAPPED>>",
                             ned_status="?", hormone_status="?", path_note=""))
        else:
            mapping[ident] = (ds, d)
            rows.append(dict(orig_ident=ident, dataset=ds, n=n, patient_phenotype=d["phenotype"],
                             ned_status=d["ned"], hormone_status=d["hormone"], path_note=d["note"]))
    tbl = pd.DataFrame(rows).sort_values(["dataset", "orig_ident"])
    print("\n[map] orig.ident -> phenotype:")
    print(tbl.to_string(index=False))

    if unmapped:
        raise SystemExit(
            f"\n[abort] {len(unmapped)} orig.ident value(s) could not be mapped: {unmapped}\n"
            f"        Check the Gao_<N>/Li_<N> numbering vs the paper patient numbers, then "
            f"edit GAO/LI dicts in this script.")

    # apply
    def col(ident, field):
        ds, d = mapping[str(ident)]
        return ds if field == "dataset" else d[field]

    if args.dataset_key not in adata.obs:
        adata.obs[args.dataset_key] = pd.Categorical(idents.map(lambda i: mapping[i][0]))
    else:
        # sanity-check existing dataset column against our derivation
        derived = idents.map(lambda i: mapping[i][0])
        mismatch = (adata.obs[args.dataset_key].astype(str).values != derived.values).sum()
        if mismatch:
            print(f"[warn] existing obs[{args.dataset_key!r}] disagrees with derived dataset "
                  f"on {mismatch} cells — keeping existing column")
    for field in ["phenotype", "ned", "hormone", "note"]:
        outcol = {"phenotype": "patient_phenotype", "ned": "ned_status",
                  "hormone": "hormone_status", "note": "path_note"}[field]
        adata.obs[outcol] = pd.Categorical(idents.map(lambda i: mapping[i][1][field]))

    # QC crosstabs
    print("\n[qc] dataset x ned_status:")
    print(pd.crosstab(adata.obs[args.dataset_key], adata.obs["ned_status"]))
    print("\n[qc] dataset x hormone_status:")
    print(pd.crosstab(adata.obs[args.dataset_key], adata.obs["hormone_status"]))

    # validation table — runs in dry-run too so you can inspect before writing
    if args.report:
        report_validation(adata, args.dataset_key, args.report_dir, use_raw=args.use_raw)

    if args.dry_run:
        print("\n[dry-run] no file written.")
        return
    out = args.out or args.adata.replace(".h5ad", ".annot.h5ad")
    adata.write_h5ad(out)
    print(f"\n[ok] wrote {out}  (+ cols: {args.dataset_key}, patient_phenotype, "
          f"ned_status, hormone_status, path_note)")


if __name__ == "__main__":
    main()
