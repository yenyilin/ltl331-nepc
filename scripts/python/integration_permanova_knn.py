#!/usr/bin/env python3
"""
integration_permanova_knn.py — global + directed integration QC for the LTL331 + Gao + Li
joint object, complementing the *local* LISI metrics in integration_lisi_kbet.py with no
scikit-bio / vegan dependency.

Two analyses
------------
(1) PERMANOVA (one-way, Anderson 2001), dataset as the factor, on the naive PCA embedding
    vs the Harmony-corrected embedding. R^2 = fraction of embedding variance explained by
    the dataset label; a DROP from PCA->Harmony quantifies how much batch structure Harmony
    removed *globally* (the complement to iLISI, which only sees local kNN mixing).
    Subsampled (distance matrix is O(n^2)); averaged over several independent subsamples.
    Also reports per-group multivariate DISPERSION (mean distance to group centroid, the
    betadisper quantity for Euclidean embeddings).

   PERMANOVA R^2 can be high for two reasons the test cannot tell apart: (1) the groups are
   centered in DIFFERENT places (real separation), or (2) they overlap in the middle but one
   is a TIGHT blob and the other a WIDE cloud (a pure spread difference). Picture two dart
   patterns: different bullseyes, versus the same bullseye but one tightly grouped and one
   scattered. A PDX (one clonal tumor, homogeneous) and clinical biopsies (many patients, heterogeneous)
   almost certainly differ in spread, so a skeptic could wave away the R^2 as "just" that.
   Reporting dispersion separately lets us show the R^2 reflects LOCATION, not just spread.

   How to read it: similar dispersions across groups + high R^2 => real separation; very
   different dispersions => interpret R^2 cautiously (it may be partly a spread artifact;
   that's what the PERMDISP/betadisper number is for).

(2) Directed kNN retrieval for the cluster-17 bridge cells: of each LTL331 cluster-17 cell's
    k nearest neighbors, what fraction are clinical (Gao/Li)? Computed on PCA (before) AND
    Harmony (after). Even when global iLISI~1, this asks the precise question we care about —
    do the bridge cells acquire clinical neighbors after integration? Reportable either way.

Directionality
  PERMANOVA R^2 : LOWER after Harmony = more batch removed (good). p tests R^2 != 0.
  kNN frac_clinical : HIGHER after Harmony = bridge cells gained clinical neighbors.

Usage
  python integration_permanova_knn.py --adata harmony.h5ad --out data/integration_qc
  python integration_permanova_knn.py --adata harmony.h5ad --query-cluster 17 --knn-k 30 \
      --permanova-n 2000 --permanova-reps 5 --permanova-perms 999
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

# default EMT signature (student's set; aligns with cl17 = EMT-mesenchymal, not basal)
EMT_MARKERS = ["CDH2", "CDH11", "FN1", "VIM", "TWIST1", "SNAI1", "ZEB1", "ZEB2", "DCN"]


# ---------------------------------------------------------------- PERMANOVA
def _ssw(D2, codes, a, ng):
    """Within-group sum of squares for a one-way grouping (vectorized, no fancy-index copy).
    SS_W = sum_g (1/n_g) * sum_{i<j in g} d^2 ; with 1_g^T D2 1_g = 2 * sum_{i<j in g} d^2."""
    G = np.zeros((codes.size, a))
    G[np.arange(codes.size), codes] = 1.0
    within = np.einsum("ij,ij->j", G, D2 @ G)        # 1_g^T D2 1_g per group
    return np.sum(0.5 * within / ng)


def permanova_oneway(X, labels, n_perm, rng):
    """One-way PERMANOVA on squared-Euclidean distances. Returns (R2, pseudo_F, p)."""
    cats, codes = np.unique(labels, return_inverse=True)
    a, N = cats.size, codes.size
    ng = np.bincount(codes, minlength=a).astype(float)
    D2 = pairwise_distances(X, metric="sqeuclidean")
    ss_t = 0.5 * D2.sum() / N
    ssw_obs = _ssw(D2, codes, a, ng)
    ssa_obs = ss_t - ssw_obs
    f_obs = (ssa_obs / (a - 1)) / (ssw_obs / (N - a))
    r2 = ssa_obs / ss_t
    ge = 1                                            # +1 for the observed (standard)
    for _ in range(n_perm):
        perm = rng.permutation(codes)
        ng_p = np.bincount(perm, minlength=a).astype(float)
        ssw_p = _ssw(D2, perm, a, ng_p)
        f_p = ((ss_t - ssw_p) / (a - 1)) / (ssw_p / (N - a))
        if f_p >= f_obs:
            ge += 1
    return float(r2), float(f_obs), (ge / (n_perm + 1))


def group_dispersion(X, labels):
    """Mean Euclidean distance to the group centroid, per group (betadisper for Euclidean X)."""
    out = {}
    for g in np.unique(labels):
        Xg = X[labels == g]
        d = np.linalg.norm(Xg - Xg.mean(axis=0), axis=1)
        out[g] = (float(d.mean()), int(Xg.shape[0]))
    return out


# ---------------------------------------------------------------- directed kNN
def knn_clinical_fraction(X, datasets, query_mask, clinical, k, emt=None):
    """For each query cell, fraction of its k NN in clinical datasets (+ per-cohort, +EMT)."""
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)        # +1 because first NN is self
    _, idx = nn.kneighbors(X[query_mask])
    idx = idx[:, 1:]                                        # drop self column
    neigh_ds = datasets[idx]                               # (n_query, k)
    is_clin = np.isin(neigh_ds, list(clinical))
    rows = {
        "n_query": int(query_mask.sum()),
        "mean_frac_clinical": float(is_clin.mean(axis=1).mean()),
        "median_frac_clinical": float(np.median(is_clin.mean(axis=1))),
    }
    for c in sorted(set(datasets)):
        rows[f"mean_frac_{c}"] = float((neigh_ds == c).mean(axis=1).mean())
    if emt is not None:
        rows["mean_neighbor_EMT"] = float(emt[idx].mean(axis=1).mean())
        rows["mean_query_EMT"] = float(emt[query_mask].mean())
    return rows


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--adata", required=True)
    p.add_argument("--out", default="data/integration_qc")
    p.add_argument("--batch-key", default="dataset")
    p.add_argument("--reps", default="X_pca,X_harmony_dataset",
                   help="comma list: naive embedding first, corrected second")
    p.add_argument("--n-dims", type=int, default=30, help="first N dims of each embedding")
    p.add_argument("--permanova-n", type=int, default=2000, help="cells per subsample")
    p.add_argument("--permanova-reps", type=int, default=5, help="independent subsamples")
    p.add_argument("--permanova-perms", type=int, default=999)
    p.add_argument("--seed", type=int, default=0)
    # directed kNN
    p.add_argument("--query-cluster-key", default="original_cluster")
    p.add_argument("--query-cluster", default="17")
    p.add_argument("--ltl331-label", default="LTL331")
    p.add_argument("--clinical-labels", default="Gao,Li")
    p.add_argument("--knn-k", type=int, default=30)
    p.add_argument("--emt-markers", default=",".join(EMT_MARKERS))
    p.add_argument("--use-raw", action="store_true", default=True)
    p.add_argument("--no-raw", dest="use_raw", action="store_false")
    p.add_argument("--no-knn", action="store_true", help="skip the directed-kNN analysis")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f"[load] {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"       {adata.n_obs:,} cells; batch {adata.obs[args.batch_key].value_counts().to_dict()}")
    reps = [r.strip() for r in args.reps.split(",") if r.strip()]
    for r in reps:
        if r not in adata.obsm:
            raise SystemExit(f"[abort] obsm[{r!r}] missing; have {list(adata.obsm)}")
    labels = adata.obs[args.batch_key].astype(str).to_numpy()

    def emb(rep):
        X = np.asarray(adata.obsm[rep])
        return X[:, : min(args.n_dims, X.shape[1])]

    # ---- (1) PERMANOVA + dispersion, before vs after --------------------
    print(f"\n[permanova] dataset effect, {args.permanova_reps}x subsample of "
          f"{args.permanova_n}, {args.permanova_perms} perms")
    pm_rows, disp_rows = [], []
    for rep in reps:
        X = emb(rep)
        r2s, fs, ps = [], [], []
        for s in range(args.permanova_reps):
            sub = rng.choice(adata.n_obs, size=min(args.permanova_n, adata.n_obs), replace=False)
            r2, f, pv = permanova_oneway(X[sub], labels[sub], args.permanova_perms, rng)
            r2s.append(r2); fs.append(f); ps.append(pv)
        pm_rows.append({"embedding": rep, "R2_mean": round(np.mean(r2s), 4),
                        "R2_sd": round(np.std(r2s), 4), "pseudoF_mean": round(np.mean(fs), 2),
                        "p_max": max(ps), "ndim": X.shape[1]})
        print(f"   {rep:22s} R2={np.mean(r2s):.4f}+-{np.std(r2s):.4f}  "
              f"F={np.mean(fs):.1f}  p<={max(ps):.3g}")
        for g, (md, n) in group_dispersion(X, labels).items():
            disp_rows.append({"embedding": rep, "group": g,
                              "mean_dist_to_centroid": round(md, 3), "n": n})
    pm = pd.DataFrame(pm_rows)
    pm.to_csv(f"{args.out}/permanova_before_after.tsv", sep="\t", index=False)
    disp = pd.DataFrame(disp_rows)
    disp.to_csv(f"{args.out}/dispersion_by_group.tsv", sep="\t", index=False)
    print("\n[dispersion] mean distance to group centroid (betadisper-style):")
    print(disp.pivot(index="group", columns="embedding",
                     values="mean_dist_to_centroid").to_string())

    delta_r2 = None
    if len(pm_rows) >= 2:
        delta_r2 = pm_rows[1]["R2_mean"] - pm_rows[0]["R2_mean"]
        print(f"\n[permanova] dR2 (corrected - naive) = {delta_r2:+.4f}  "
              f"({'batch structure REDUCED' if delta_r2 < 0 else 'NOT reduced'})")

    # ---- (2) directed kNN retrieval for cluster-17 ----------------------
    knn_rows = []
    if not args.no_knn:
        clinical = set(args.clinical_labels.split(","))
        qkey = args.query_cluster_key
        if qkey not in adata.obs:
            print(f"[knn] WARN obs[{qkey!r}] absent ({list(adata.obs)[:8]}...) — skipping kNN")
        else:
            qmask = ((labels == args.ltl331_label)
                     & (adata.obs[qkey].astype(str).to_numpy() == str(args.query_cluster)))
            if qmask.sum() == 0:
                print(f"[knn] WARN no {args.ltl331_label} cells with {qkey}=={args.query_cluster}")
            else:
                emt = None
                genes = [g for g in args.emt_markers.split(",") if g]
                src = (adata.raw.var_names if (args.use_raw and adata.raw is not None)
                       else adata.var_names)
                present = [g for g in genes if g in set(src)]
                if present:
                    sc.tl.score_genes(adata, present, score_name="_emt", use_raw=args.use_raw)
                    emt = adata.obs["_emt"].to_numpy()
                    print(f"\n[knn] EMT signature {len(present)}/{len(genes)} genes present")
                print(f"[knn] {int(qmask.sum())} query cells "
                      f"({args.ltl331_label} {qkey}={args.query_cluster}), k={args.knn_k}")
                for rep in reps:
                    row = {"embedding": rep,
                           **knn_clinical_fraction(emb(rep), labels, qmask, clinical,
                                                    args.knn_k, emt)}
                    knn_rows.append(row)
                    print(f"   {rep:22s} frac_clinical(mean)={row['mean_frac_clinical']:.3f} "
                          f"(Gao {row.get('mean_frac_Gao', 0):.3f}, Li {row.get('mean_frac_Li', 0):.3f})")
                pd.DataFrame(knn_rows).to_csv(
                    f"{args.out}/cluster{args.query_cluster}_knn_retrieval.tsv",
                    sep="\t", index=False)

    # ---- verdict --------------------------------------------------------
    verdict = {"permanova": pm_rows, "delta_R2_corrected_minus_naive": delta_r2,
               "dispersion": disp_rows, "knn_retrieval": knn_rows,
               "note": "R2 lower after Harmony = global batch reduced; iLISI (see "
                       "lisi_kbet_metrics.tsv) shows local kNN mixing separately."}
    with open(f"{args.out}/permanova_knn_verdict.json", "w") as fh:
        json.dump(verdict, fh, indent=2)
    print(f"\n[ok] wrote {args.out}/permanova_before_after.tsv, dispersion_by_group.tsv, "
          f"cluster{args.query_cluster}_knn_retrieval.tsv, permanova_knn_verdict.json")


if __name__ == "__main__":
    main()
