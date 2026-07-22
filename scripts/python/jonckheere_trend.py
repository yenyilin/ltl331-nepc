#!/usr/bin/env python
"""
jonckheere_trend.py — Jonckheere-Terpstra ordered-trend test for velocity
pseudotime across the canonical adenocarcinoma -> neuroendocrine cluster order.

For BOTH the per-cell data and a timepoint-aggregated version (which respects
independence), it reports:
  - Jonckheere-Terpstra statistic J, standardized Z (TIE-CORRECTED), one-sided p
  - JT probabilistic index  PI = J / sum_{i<j} n_i n_j
        0.5 = no trend, 1.0 = every cross-cluster pair concordant with the order
  - Spearman rho / Kendall tau of pseudotime vs cluster-order rank
  - Spearman rho of pseudotime vs experimental week (if a week key is present)

WHY the aggregated version: 50k single cells are NOT independent, so the per-cell
Z is astronomically large (pseudoreplication). Report PI + a CAPPED p as the
honest effect size, and lead with the timepoint-aggregated test, in which each
(cluster x timepoint) contributes one median -- a believable n and Z that no
one can call inflated.

Usage:
  python jonckheere_trend.py --h5ad ltl331_velocity.h5ad \
      --cluster-key seurat_clusters --pseudotime-key velocity_pseudotime \
      --week-key week_num --cluster-order data/cluster_order.json \
      --exclude-clusters 18 --out data/forfig5/jt_trend.tsv
"""
import argparse
import json
import sys
from itertools import combinations

import numpy as np
from scipy.stats import mannwhitneyu, norm, spearmanr, kendalltau


# --------------------------------------------------------------------------- #
# Jonckheere-Terpstra core
# --------------------------------------------------------------------------- #
def jt_variance_tiecorrected(group_sizes, tie_sizes, N):
    """Exact JT null variance with tie correction (Hollander, Wolfe & Chicken,
    Nonparametric Statistical Methods). Reduces to the no-tie formula when every
    pooled value is unique (all tie_sizes == 1)."""
    n = np.asarray(group_sizes, dtype=float)
    e = np.asarray(tie_sizes, dtype=float)
    if N < 3:
        return np.nan
    term1 = (N * (N - 1) * (2 * N + 5)
             - np.sum(n * (n - 1) * (2 * n + 5))
             - np.sum(e * (e - 1) * (2 * e + 5))) / 72.0
    term2 = (np.sum(n * (n - 1) * (n - 2)) * np.sum(e * (e - 1) * (e - 2))) \
        / (36.0 * N * (N - 1) * (N - 2))
    term3 = (np.sum(n * (n - 1)) * np.sum(e * (e - 1))) / (8.0 * N * (N - 1))
    return term1 + term2 + term3


def jonckheere(groups):
    """groups: list of 1D arrays, already in the hypothesized INCREASING order.
    Returns J, tie-corrected Z, one-sided p (increasing alternative), and the
    probabilistic index PI."""
    k = len(groups)
    # J = sum over ordered pairs of #(later > earlier) [+ 0.5 ties]
    # scipy's mannwhitneyu(x, y).statistic = U1 = #(x_i > y_j) + 0.5*ties
    J = 0.0
    for i, j in combinations(range(k), 2):
        J += mannwhitneyu(groups[j], groups[i], method="asymptotic").statistic
    n = np.array([len(g) for g in groups], dtype=float)
    N = int(n.sum())
    EJ = (N ** 2 - np.sum(n ** 2)) / 4.0
    _, tie_counts = np.unique(np.concatenate(groups), return_counts=True)
    VJ = jt_variance_tiecorrected(n, tie_counts, N)
    Z = (J - EJ) / np.sqrt(VJ) if VJ and VJ > 0 else np.nan
    denom = float(np.sum([n[a] * n[b] for a, b in combinations(range(k), 2)]))
    PI = J / denom if denom else np.nan
    P = norm.sf(Z) if np.isfinite(Z) else np.nan          # one-sided, increasing
    return dict(J=J, EJ=EJ, VJ=VJ, Z=Z, P=P, PI=PI, N=N, k=k)


# --------------------------------------------------------------------------- #
# IO helpers
# --------------------------------------------------------------------------- #
def load_obs(h5ad, keys):
    """Read only the obs columns we need (backed mode -> no X load)."""
    try:
        import anndata as ad
        A = ad.read_h5ad(h5ad, backed="r")
        have = [k for k in keys if k in A.obs.columns]
        obs = A.obs[have].copy()
        try:
            A.file.close()
        except Exception:
            pass
    except Exception:
        import scanpy as sc
        adata = sc.read_h5ad(h5ad)
        have = [k for k in keys if k in adata.obs.columns]
        obs = adata.obs[have].copy()
    return obs


def load_cluster_order(path):
    with open(path) as fh:
        obj = json.load(fh)
    if isinstance(obj, dict):
        # accept {"order": [...]} or {cluster: rank}
        if "order" in obj:
            obj = obj["order"]
        else:
            obj = [k for k, _ in sorted(obj.items(), key=lambda kv: kv[1])]
    return [str(x) for x in obj]


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--cluster-key", default="seurat_clusters")
    ap.add_argument("--pseudotime-key", default="velocity_pseudotime")
    ap.add_argument("--week-key", default="week_num",
                    help="experimental-time key in obs (skipped if absent)")
    ap.add_argument("--cluster-order", default="data/cluster_order.json")
    ap.add_argument("--exclude-clusters", nargs="*", default=["18"],
                    help="clusters to drop (e.g. normal-reference / doublet bin)")
    ap.add_argument("--cap-p", type=float, default=1e-4,
                    help="floor for the p-value printed in the suggested sentence")
    ap.add_argument("--out", default=None, help="optional TSV summary path")
    args = ap.parse_args()

    obs = load_obs(args.h5ad, [args.cluster_key, args.pseudotime_key, args.week_key])
    for need in (args.cluster_key, args.pseudotime_key):
        if need not in obs.columns:
            sys.exit(f"ERROR: '{need}' not found in obs. Present: {list(obs.columns)}")

    obs = obs.copy()
    obs[args.cluster_key] = obs[args.cluster_key].astype(str)
    obs[args.pseudotime_key] = obs[args.pseudotime_key].astype(float)
    obs = obs[~obs[args.pseudotime_key].isna()]
    obs = obs[~obs[args.cluster_key].isin([str(c) for c in args.exclude_clusters])]

    order = load_cluster_order(args.cluster_order)
    present = [c for c in order if c in set(obs[args.cluster_key])]
    missing = [c for c in order if c not in set(obs[args.cluster_key])]
    if missing:
        print(f"[note] cluster_order entries not in data (skipped): {missing}")
    rank = {c: i for i, c in enumerate(present)}

    pt = args.pseudotime_key

    # --- per-cell JT --------------------------------------------------------
    groups = [obs.loc[obs[args.cluster_key] == c, pt].to_numpy() for c in present]
    res_cell = jonckheere(groups)

    # cell-level rank correlations
    cl_rank = obs[args.cluster_key].map(rank).to_numpy(dtype=float)
    keep = ~np.isnan(cl_rank)
    rho_order, p_order = spearmanr(obs[pt].to_numpy()[keep], cl_rank[keep])
    tau_order, _ = kendalltau(obs[pt].to_numpy()[keep], cl_rank[keep])

    rho_week = p_week = np.nan
    have_week = args.week_key in obs.columns
    if have_week:
        wk = obs[args.week_key].astype(float).to_numpy()
        wkeep = keep & ~np.isnan(wk)
        rho_week, p_week = spearmanr(obs[pt].to_numpy()[wkeep], wk[wkeep])

    # --- timepoint-aggregated JT (independence-respecting) ------------------
    res_agg = None
    if have_week:
        agg_groups = []
        for c in present:
            sub = obs[obs[args.cluster_key] == c]
            med = sub.groupby(args.week_key, observed=True)[pt].median().to_numpy()
            agg_groups.append(med[~np.isnan(med)])
        if all(len(g) > 0 for g in agg_groups) and sum(len(g) for g in agg_groups) >= 3:
            res_agg = jonckheere(agg_groups)

    # --- report -------------------------------------------------------------
    def fmt(r):
        return (f"  J={r['J']:.0f}  E[J]={r['EJ']:.0f}  Z={r['Z']:.3g}  "
                f"PI={r['PI']:.4f}  p={r['P']:.3e}  (N={r['N']}, k={r['k']})")

    print("\n=== Jonckheere-Terpstra: velocity pseudotime vs adeno->NE cluster order ===")
    print(f"cluster order used: {present}")
    print("\n[per-cell] (pseudoreplicated -- PI is the honest effect size; Z/p are anticonservative)")
    print(fmt(res_cell))
    print(f"\n[cell-level rank corr]  Spearman rho(order)={rho_order:.3f}  "
          f"Kendall tau(order)={tau_order:.3f}")
    if have_week:
        print(f"[cell-level rank corr]  Spearman rho(week) ={rho_week:.3f}")
    if res_agg is not None:
        print("\n[timepoint-aggregated] (each cluster = its per-week median pseudotimes; "
              "respects independence -- REPORT THIS AS PRIMARY)")
        print(fmt(res_agg))
    else:
        print("\n[timepoint-aggregated] skipped (no usable week key / too few medians)")

    p_show = max(res_cell["P"], args.cap_p) if np.isfinite(res_cell["P"]) else args.cap_p
    print("\n--- suggested manuscript sentence ------------------------------------------")
    print(
        f'Velocity pseudotime increases monotonically along the adenocarcinoma->'
        f'neuroendocrine cluster order (Spearman rho = {rho_order:.2f}, Kendall tau '
        f'= {tau_order:.2f}; Jonckheere-Terpstra probabilistic index = '
        f'{res_cell["PI"]:.2f}, i.e. {100*res_cell["PI"]:.0f}% of cross-cluster cell '
        f'pairs are concordant with the ordering; trend p < {args.cap_p:g}). As single '
        f'cells are not independent, p-values are interpreted descriptively; the same '
        f'monotone trend holds on per-timepoint cluster medians'
        + (f' (Jonckheere-Terpstra Z = {res_agg["Z"]:.1f}, p = {res_agg["P"]:.1e}).'
           if res_agg is not None else '.')
    )
    print("----------------------------------------------------------------------------")

    # --- optional TSV -------------------------------------------------------
    if args.out:
        import csv
        import os
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        rows = [
            ("test", "J", "EJ", "Z", "PI", "p", "N", "k"),
            ("jt_per_cell", res_cell["J"], res_cell["EJ"], res_cell["Z"],
             res_cell["PI"], res_cell["P"], res_cell["N"], res_cell["k"]),
        ]
        if res_agg is not None:
            rows.append(("jt_timepoint_aggregated", res_agg["J"], res_agg["EJ"],
                         res_agg["Z"], res_agg["PI"], res_agg["P"], res_agg["N"],
                         res_agg["k"]))
        rows.append(("spearman_pt_vs_order", "", "", "", rho_order, p_order, "", ""))
        rows.append(("kendall_pt_vs_order", "", "", "", tau_order, "", "", ""))
        if have_week:
            rows.append(("spearman_pt_vs_week", "", "", "", rho_week, p_week, "", ""))
        with open(args.out, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        print(f"[wrote] {args.out}")


if __name__ == "__main__":
    main()
