#!/usr/bin/env python3
"""
plot_fig6a_label_transfer.py — Fig 6A: project clinical cells onto LTL331 coarse states by
reference projection (NO co-embedding), then summarize at the stable patient level and show a
shared phenotype-axis co-map. Companion to clinical_signature_validation.py (which gives the
continuous patient x signature heatmap); this adds the discrete label-transfer + phenotype map.

WHY reference projection and NOT Harmony/scANVI label transfer:
  Joint co-embedding fails here (iLISI~1.0 — datasets do not mix), which also breaks anchor-based
  transfer that assumes a shared manifold. Reference projection (centroid correlation / signature
  argmax) is the method class built for the no-co-embedding case. Two independent projectors are
  run and their AGREEMENT is reported as robustness.

WHY patient-level and coarse classes:
  Per-dataset clinical cluster IDs (Gao_0/Gao_20/Li_5/...) move every time clustering is redone, so
  they are NOT a stable validation unit. Patient identity + biopsy pathology are stable. Coarse
  4-class labels (PRAD / Intermediate / NE-ASCL1+ / NE-ASCL1-) transfer far more reliably across
  cohorts than 18 fine clusters. Use --merge-ne for the conservative 3-class (PRAD/Intermediate/NEPC)
  view; the ASCL1+/- split in clinical biopsies is EXPLORATORY (state so in the legend).

Outputs (in --out):
  fig6a_coarse_class_signatures.tsv    reference-derived marker genes per coarse class
  fig6a_cell_labels.tsv                per clinical cell: both projections + agreement flag
  fig6a_patient_composition.tsv        per patient: transferred-class fractions + pathology
  fig6a_concordance_stats.tsv          NE-fraction vs NE pathology (MWU + AUC), LOPO + subsample
  fig6a_composition.{pdf,png}          patient composition stacked bar (disease-axis ordered)
  fig6a_phenotype_map.{pdf,png}        PRAD-score x NE-score; LTL331 density + clinical overlay

Pure scanpy/numpy/pandas/scipy/matplotlib. Reads lognorm from .raw (LTL331 .X is scaled).

Usage
  python plot_fig6a_label_transfer.py --h5ad harmony.annotated.h5ad \\
      --group-col dataset --ltl331-label LTL331 --patient-col orig.ident \\
      --phenotype-col patient_phenotype --cluster-col seurat_clusters \\
      --cluster-order data/cluster_order.json --out data/fig6a --plot-format pdf png
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig  # noqa: E402

# Default coarse class -> LTL331 seurat_clusters. EDIT/override via --class-map JSON.
# Derived from cluster_order.json stages; ASCL1 split: cl10 ASCL1+, cl7/13/1/3/9/14 ASCL1-.
DEFAULT_CLASS_MAP = {
    "PRAD":       ["11", "15", "0", "2", "5", "8", "12", "16", "4", "6"],
    "Intermediate": ["17"],
    "NE-ASCL1+":  ["10"],
    "NE-ASCL1-":  ["13", "1", "3", "9", "14", "7"],
}
CLASS_ORDER = ["PRAD", "Intermediate", "NE-ASCL1+", "NE-ASCL1-"]
CLASS_COLORS = {"PRAD": "#2c7fb8", "Intermediate": "#fdae61",
                "NE-ASCL1+": "#b2182b", "NE-ASCL1-": "#762a83", "NEPC": "#b2182b"}

# patient_phenotype -> (disease_state, disease_axis ordinal). Same scheme as
# clinical_signature_validation.py; override with --metadata.
DEFAULT_AXIS = {
    "Gao_PIN": ("PIN", 0), "Li_mHSPC_adeno": ("HSPC", 1),
    "Li_nmHSPC_mixedNEPC": ("HSPC", 1), "Gao_PRAD": ("PRAD", 2),
    "Li_mCRPC_mixedNEPC": ("CRPC", 3), "Gao_NEPC_smallcell": ("NEPC", 4),
    "Li_mCRPC_pureSmallCell": ("NEPC", 4),
}


def dense(x):
    return x.toarray() if ssp.issparse(x) else np.asarray(x)


def raw_frame(adata, genes):
    """lognorm expression (cells x genes) from .raw, restricted to `genes` present there."""
    src = adata.raw if adata.raw is not None else adata
    var = list(map(str, src.var_names))
    idx = {g: i for i, g in enumerate(var)}
    keep = [g for g in genes if g in idx]
    M = dense(src.X[:, [idx[g] for g in keep]])
    return pd.DataFrame(M, index=adata.obs_names, columns=keep)


def per_cell_pearson(Q, centroids):
    """Q: (n,G) z-scored rows; centroids: dict label->(G,) z-scored. Returns DataFrame n x labels."""
    out = {}
    for lab, c in centroids.items():
        out[lab] = (Q * c).sum(axis=1) / Q.shape[1]
    return pd.DataFrame(out, index=range(Q.shape[0]))


def zscore_rows(M):
    mu = M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (M - mu) / sd


def auc(scores, labels):
    """ROC AUC of `scores` predicting binary `labels` (1=positive). Rank-based, ties averaged."""
    s = np.asarray(scores, float); y = np.asarray(labels, int)
    n1 = y.sum(); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(s, kind="stable")
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    df = pd.Series(s).rank(method="average").to_numpy()
    return (df[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True, help="integrated object (LTL331 + clinical); .raw lognorm")
    ap.add_argument("--group-col", default="dataset")
    ap.add_argument("--ltl331-label", default="LTL331")
    ap.add_argument("--patient-col", default="orig.ident")
    ap.add_argument("--phenotype-col", default="patient_phenotype")
    ap.add_argument("--ned-col", default="ned_status")
    ap.add_argument("--cluster-col", default="seurat_clusters")
    ap.add_argument("--cluster-order", default=None)  # kept for parity; class map drives grouping
    ap.add_argument("--class-map", default=None, help="JSON {class: [cluster ids]} overriding default")
    ap.add_argument("--merge-ne", action="store_true",
                    help="collapse NE-ASCL1+/- into a single NEPC class (conservative 3-class)")
    ap.add_argument("--metadata", default=None, help="TSV phenotype->(disease_state,disease_axis)")
    ap.add_argument("--top-n", type=int, default=40, help="marker genes per class for signatures")
    ap.add_argument("--n-subsample", type=int, default=50, help="cell-subsample reps for robustness")
    ap.add_argument("--subsample-frac", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="data/fig6a")
    ap.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f"[load] {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    grp = adata.obs[args.group_col].astype(str)
    is_ref = grp == args.ltl331_label
    ref = adata[is_ref].copy()
    qry = adata[~is_ref].copy()
    print(f"  reference (LTL331): {ref.n_obs}   clinical query: {qry.n_obs}")

    # ---- coarse class on reference ----
    class_map = json.load(open(args.class_map)) if args.class_map else dict(DEFAULT_CLASS_MAP)
    classes = list(CLASS_ORDER)
    if args.merge_ne:
        class_map = {"PRAD": class_map["PRAD"], "Intermediate": class_map["Intermediate"],
                     "NEPC": class_map["NE-ASCL1+"] + class_map["NE-ASCL1-"]}
        classes = ["PRAD", "Intermediate", "NEPC"]
    cl2class = {str(c): k for k, cs in class_map.items() for c in cs}
    ref_cl = ref.obs[args.cluster_col].astype(str)
    ref.obs["coarse"] = ref_cl.map(cl2class)
    print("[class map] (EDIT via --class-map if ASCL1 split differs):")
    for k in classes:
        print(f"   {k:14s} <- clusters {class_map[k]}  ({(ref.obs['coarse']==k).sum()} cells)")
    ref = ref[ref.obs["coarse"].notna()].copy()
    ref.obs["coarse"] = pd.Categorical(ref.obs["coarse"], categories=classes, ordered=True)

    # ---- reference-derived class signatures (rank_genes_groups on .raw lognorm) ----
    use_raw = ref.raw is not None
    sc.tl.rank_genes_groups(ref, "coarse", method="wilcoxon", use_raw=use_raw, n_genes=args.top_n)
    sig = {}
    for k in classes:
        sig[k] = [g for g in sc.get.rank_genes_groups_df(ref, group=k)["names"].tolist()[:args.top_n]]
    pd.DataFrame({k: pd.Series(v) for k, v in sig.items()}).to_csv(
        f"{args.out}/fig6a_coarse_class_signatures.tsv", sep="\t", index=False)
    feat = sorted({g for v in sig.values() for g in v})

    # restrict feature genes to those present in BOTH ref.raw and qry.raw
    rsrc = set(map(str, (ref.raw.var_names if use_raw else ref.var_names)))
    qsrc = set(map(str, (qry.raw.var_names if qry.raw is not None else qry.var_names)))
    feat = [g for g in feat if g in rsrc and g in qsrc]
    print(f"[features] {len(feat)} marker genes shared by reference and clinical")

    # ---- projector 1: centroid correlation over the marker feature space ----
    Rf = raw_frame(ref, feat)[feat]
    centroid_raw = {k: Rf.loc[ref.obs["coarse"] == k].mean(axis=0).to_numpy() for k in classes}
    centroid = {k: (v - v.mean()) / (v.std() if v.std() else 1.0) for k, v in centroid_raw.items()}
    Qf = raw_frame(qry, feat)[feat].to_numpy()
    Qz = zscore_rows(Qf)
    corr = per_cell_pearson(Qz, centroid)
    proj_corr = corr.idxmax(axis=1).to_numpy()

    # ---- projector 2: signature argmax (sc.tl.score_genes per class) ----
    score_cols = []
    for k in classes:
        present = [g for g in sig[k] if g in qsrc]
        sc.tl.score_genes(qry, present, score_name=f"score_{k}", use_raw=(qry.raw is not None))
        score_cols.append(f"score_{k}")
    S = qry.obs[score_cols].copy()
    proj_sig = S.to_numpy().argmax(axis=1)
    proj_sig = np.array([classes[i] for i in proj_sig])

    agree = proj_corr == proj_sig
    print(f"[agreement] centroid vs signature argmax: {agree.mean():.1%} of clinical cells")

    labels = pd.DataFrame({
        args.patient_col: qry.obs[args.patient_col].astype(str).to_numpy(),
        "proj_centroid": proj_corr, "proj_signature": proj_sig, "agree": agree,
        "corr_max": corr.max(axis=1).to_numpy(),
    }, index=qry.obs_names)
    for c in score_cols:
        labels[c] = qry.obs[c].to_numpy()
    labels.to_csv(f"{args.out}/fig6a_cell_labels.tsv", sep="\t")

    # ---- patient-level composition (primary projector = centroid) ----
    axis_map = dict(DEFAULT_AXIS)
    if args.metadata and os.path.exists(args.metadata):
        m = pd.read_csv(args.metadata, sep=None, engine="python")
        axis_map.update({r[args.phenotype_col]: (r["disease_state"], int(r["disease_axis"]))
                         for _, r in m.iterrows()})
    pobs = qry.obs.copy()
    pobs["proj"] = proj_corr
    comp = (pobs.groupby(args.patient_col, observed=True)["proj"]
            .value_counts(normalize=True).unstack(fill_value=0.0).reindex(columns=classes, fill_value=0.0))
    phen = pobs.groupby(args.patient_col, observed=True)[args.phenotype_col].agg(
        lambda s: s.astype(str).mode().iat[0])
    comp["phenotype"] = phen
    comp["disease_state"] = comp["phenotype"].map(lambda p: axis_map.get(p, ("?", 99))[0])
    comp["disease_axis"] = comp["phenotype"].map(lambda p: axis_map.get(p, ("?", 99))[1])
    if args.ned_col in pobs:
        comp["ned_status"] = pobs.groupby(args.patient_col, observed=True)[args.ned_col].agg(
            lambda s: s.astype(str).mode().iat[0])
    comp = comp.sort_values(["disease_axis", "disease_state"])
    comp.to_csv(f"{args.out}/fig6a_patient_composition.tsv", sep="\t")

    # ---- concordance: per-patient NE fraction predicts NE pathology ----
    ne_classes = ["NEPC"] if args.merge_ne else ["NE-ASCL1+", "NE-ASCL1-"]
    comp["NE_fraction"] = comp[ne_classes].sum(axis=1)
    if "ned_status" in comp:
        ne_path = (comp["ned_status"].astype(str) == "NED").astype(int)
    else:
        ne_path = (comp["disease_state"] == "NEPC").astype(int)
    from scipy.stats import mannwhitneyu
    a = comp.loc[ne_path == 1, "NE_fraction"]; b = comp.loc[ne_path == 0, "NE_fraction"]
    full_auc = auc(comp["NE_fraction"], ne_path)
    mwu_p = mannwhitneyu(a, b, alternative="greater")[1] if len(a) and len(b) else np.nan

    # leave-one-patient-out AUC stability
    lopo = []
    for p in comp.index:
        sub = comp.drop(index=p)
        lopo.append(auc(sub["NE_fraction"], ne_path.drop(index=p)))
    lopo = np.array([x for x in lopo if not np.isnan(x)])
    # cell-subsample stability (recompute composition AUC on a fraction of clinical cells)
    sub_aucs = []
    n = qry.n_obs
    for _ in range(args.n_subsample):
        sel = rng.random(n) < args.subsample_frac
        sp = pobs.loc[sel]
        c2 = (sp.groupby(args.patient_col, observed=True)["proj"].value_counts(normalize=True)
              .unstack(fill_value=0.0).reindex(columns=classes, fill_value=0.0))
        nef = c2[ne_classes].sum(axis=1).reindex(comp.index).fillna(0.0)
        sub_aucs.append(auc(nef, ne_path))
    sub_aucs = np.array([x for x in sub_aucs if not np.isnan(x)])
    stats = pd.DataFrame([{
        "projector_agreement": round(float(agree.mean()), 4),
        "NEfrac_AUC_full": round(float(full_auc), 4),
        "NEfrac_MWU_p_greater": f"{mwu_p:.3g}",
        "LOPO_AUC_min": round(float(lopo.min()), 4) if len(lopo) else np.nan,
        "LOPO_AUC_median": round(float(np.median(lopo)), 4) if len(lopo) else np.nan,
        "subsample_AUC_median": round(float(np.median(sub_aucs)), 4) if len(sub_aucs) else np.nan,
        "subsample_AUC_p05": round(float(np.quantile(sub_aucs, 0.05)), 4) if len(sub_aucs) else np.nan,
    }])
    stats.to_csv(f"{args.out}/fig6a_concordance_stats.tsv", sep="\t", index=False)
    print("[concordance]\n" + stats.to_string(index=False))

    # ---- plot 1: patient composition stacked bar ----
    plt = setup_style(font_size=7)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.42 * len(comp) + 3.0), 4.2))
    bottom = np.zeros(len(comp))
    x = np.arange(len(comp))
    for k in classes:
        ax.bar(x, comp[k].to_numpy(), bottom=bottom, width=0.8,
               color=CLASS_COLORS.get(k, "#888"), label=k, edgecolor="white", lw=0.3)
        bottom += comp[k].to_numpy()
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\n[{comp.loc[p,'disease_state']}]" for p in comp.index],
                       rotation=90, fontsize=6)
    ax.set_ylabel("fraction of clinical cells (transferred LTL331 state)")
    ax.set_ylim(0, 1)
    prev = None
    for i, st in enumerate(comp["disease_state"]):
        if prev is not None and st != prev:
            ax.axvline(i - 0.5, color="#333", lw=0.7)
        prev = st
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=6, frameon=False,
              title="LTL331 state")
    ax.set_title("Clinical patients: LTL331 state composition by reference projection", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    savefig(fig, args.out, "fig6a_composition", args.plot_format)
    plt.close(fig)

    # ---- plot 2: phenotype-axis co-map (LTL331 density backdrop + clinical overlay) ----
    # axes: PRAD-class score (x) vs combined-NE score (y); honest substitute for co-embedding.
    for k in classes:
        present = [g for g in sig[k] if g in rsrc]
        sc.tl.score_genes(ref, present, score_name=f"score_{k}", use_raw=use_raw)
    ne_ref = ref.obs[[f"score_{k}" for k in ne_classes]].max(axis=1)
    ne_qry = qry.obs[[f"score_{k}" for k in ne_classes]].max(axis=1)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.hexbin(ref.obs["score_PRAD"], ne_ref, gridsize=45, cmap="Greys", bins="log",
              mincnt=1, zorder=1)
    for k in classes:
        m = proj_corr == k
        ax.scatter(qry.obs["score_PRAD"].to_numpy()[m], ne_qry.to_numpy()[m], s=6,
                   c=CLASS_COLORS.get(k, "#888"), alpha=0.6, edgecolor="none",
                   label=f"clinical → {k}", zorder=3)
    ax.set_xlabel("PRAD-state signature score")
    ax.set_ylabel("NE-state signature score (max of NE classes)")
    ax.set_title("Shared phenotype space: LTL331 (grey density) + clinical cells\n"
                 "(signature-score axes — no batch integration)", fontsize=8)
    ax.legend(loc="upper right", fontsize=6, frameon=False, markerscale=1.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    savefig(fig, args.out, "fig6a_phenotype_map", args.plot_format)
    plt.close(fig)

    print(f"[ok] wrote fig6a_* (labels, composition, concordance, 2 figures) to {args.out}/")
    if args.merge_ne:
        print("[note] --merge-ne: NE shown as one class; ASCL1+/- split omitted (conservative).")
    else:
        print("[note] ASCL1+/- split in CLINICAL cells is exploratory — label it so in the caption.")


if __name__ == "__main__":
    main()
