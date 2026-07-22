#!/usr/bin/env python3
"""
doublet_unify.py — derive all per-cluster / per-timepoint doublet counts from the two
harmonized count matrices and flag any data-integrity anomalies.

Inputs (both TSV, IDENTICAL layout — timepoint rows x cluster_N columns):
  counts_original.tsv : cells per (timepoint, cluster) BEFORE DoubletFinder.
  counts_valid.tsv    : the same grid AFTER doublet removal.

These two files carry identical timepoint labels (PreCx, wk2 … wk22) and identical
cluster columns (cluster_0 … cluster_17), so no format reconciliation is needed — the
script simply differences them (doublets = original − valid). Any cell where valid >
original (which would imply cells created by removal) is reported in unify_flags.txt
rather than silently emitting a negative rate. The two matrices are authoritative; the
older per_cluster_*_doublet_rates.tsv derivations disagreed with them (their "original"
column was mis-sourced, producing spurious negatives) and are NOT used here.

Outputs (OUT dir):
  matrix_original.tsv          clean wide  timepoint x cluster  (before)
  matrix_valid.tsv             clean wide  timepoint x cluster  (after)
  doublet_counts_long.tsv      tidy: timepoint, tp_order, cluster, original_n, valid_n,
                                      doublets, doublet_pct
  doublet_by_cluster.tsv       cluster, original_n, valid_n, doublets, doublet_pct,
                                      original_frac, valid_frac, frac_delta_pp[, manuscript_role]
  doublet_by_timepoint.tsv     timepoint, original_n, valid_n, doublets, doublet_pct
  unify_flags.txt              grand totals + any anomalies (negative doublets, label mismatch)

Usage:
  python doublet_unify.py --dir data/doubletfinder --out data/doubletfinder/unified
"""
import argparse
import os

import numpy as np
import pandas as pd

TP_ORDER = ["PreCx", "wk2", "wk4", "wk8", "wk12", "wk16", "wk20", "wk22"]


def load_counts(path):
    """Load a harmonized count matrix (timepoint index, cluster_N columns) -> int cluster cols."""
    m = pd.read_csv(path, sep="\t").set_index("timepoint")
    m.columns = [int(str(c).replace("cluster_", "")) for c in m.columns]
    m.index = [str(i) for i in m.index]
    m.index.name = "timepoint"
    return m.sort_index(axis=1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default="data/doubletfinder",
                   help="dir holding counts_original.tsv + counts_valid.tsv")
    p.add_argument("--original", default=None, help="override path to counts_original.tsv")
    p.add_argument("--valid", default=None, help="override path to counts_valid.tsv")
    p.add_argument("--roles", default=None,
                   help="optional TSV with columns cluster,manuscript_role to annotate "
                        "doublet_by_cluster.tsv (off by default)")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    orig_path = args.original or os.path.join(args.dir, "counts_original.tsv")
    valid_path = args.valid or os.path.join(args.dir, "counts_valid.tsv")
    out = args.out or os.path.join(args.dir, "unified")
    os.makedirs(out, exist_ok=True)
    flags = []

    O = load_counts(orig_path)
    V = load_counts(valid_path)
    print(f"[load] original {O.shape} (tp x cluster) | valid {V.shape}")

    # ---- the two grids must be identically annotated ----------------------
    if list(O.index) != list(V.index):
        flags.append(f"timepoint labels differ: original={list(O.index)} valid={list(V.index)}")
    if list(O.columns) != list(V.columns):
        flags.append(f"cluster columns differ: original={list(O.columns)} valid={list(V.columns)}")

    # ---- order rows canonically (adeno->NE timeline), keep any extras ------
    tps = [t for t in TP_ORDER if t in O.index or t in V.index]
    extra = sorted((set(O.index) | set(V.index)) - set(tps))
    if extra:
        flags.append(f"timepoints not in canonical order, appended: {extra}")
        tps += extra
    clusters = sorted(set(O.columns) | set(V.columns))
    O = O.reindex(index=tps, columns=clusters).fillna(0).astype(int)
    V = V.reindex(index=tps, columns=clusters).fillna(0).astype(int)
    O.to_csv(f"{out}/matrix_original.tsv", sep="\t")
    V.to_csv(f"{out}/matrix_valid.tsv", sep="\t")

    # ---- tidy long --------------------------------------------------------
    long = (O.reset_index().melt("timepoint", var_name="cluster", value_name="original_n")
            .merge(V.reset_index().melt("timepoint", var_name="cluster", value_name="valid_n"),
                   on=["timepoint", "cluster"]))
    long["doublets"] = long["original_n"] - long["valid_n"]
    long["doublet_pct"] = np.where(long["original_n"] > 0,
                                   (100 * long["doublets"] / long["original_n"]).round(2), 0.0)
    long["tp_order"] = long["timepoint"].map({t: i for i, t in enumerate(tps)})
    long = long.sort_values(["cluster", "tp_order"])[
        ["timepoint", "tp_order", "cluster", "original_n", "valid_n", "doublets", "doublet_pct"]]
    long.to_csv(f"{out}/doublet_counts_long.tsv", sep="\t", index=False)

    # ---- by cluster -------------------------------------------------------
    bc = pd.DataFrame({"original_n": O.sum(0), "valid_n": V.sum(0)})
    bc.index.name = "cluster"
    bc["doublets"] = bc["original_n"] - bc["valid_n"]
    bc["doublet_pct"] = np.where(bc["original_n"] > 0,
                                 (100 * bc["doublets"] / bc["original_n"]).round(2), 0.0)
    bc["original_frac"] = (bc["original_n"] / bc["original_n"].sum()).round(4)
    bc["valid_frac"] = (bc["valid_n"] / bc["valid_n"].sum()).round(4)
    bc["frac_delta_pp"] = (100 * (bc["valid_frac"] - bc["original_frac"])).round(3)
    if args.roles and os.path.exists(args.roles):
        try:
            roles = pd.read_csv(args.roles, sep="\t").set_index("cluster")["manuscript_role"]
            bc["manuscript_role"] = bc.index.map(roles)
        except Exception as e:
            flags.append(f"could not merge manuscript_role from {args.roles}: {e}")
    bc.reset_index().to_csv(f"{out}/doublet_by_cluster.tsv", sep="\t", index=False)

    # ---- by timepoint -----------------------------------------------------
    bt = pd.DataFrame({"original_n": O.sum(1), "valid_n": V.sum(1)})
    bt["doublets"] = bt["original_n"] - bt["valid_n"]
    bt["doublet_pct"] = np.where(bt["original_n"] > 0,
                                 (100 * bt["doublets"] / bt["original_n"]).round(2), 0.0)
    bt.index.name = "timepoint"
    bt.reset_index().to_csv(f"{out}/doublet_by_timepoint.tsv", sep="\t", index=False)

    # ---- integrity flags --------------------------------------------------
    neg = bc.index[bc["doublets"] < 0].tolist()
    if neg:
        flags.append(f"NEGATIVE doublets (valid>original) in clusters {neg} — "
                     f"check the input matrices; those rates were clipped to ~0 downstream.")
    grand = (int(O.values.sum()), int(V.values.sum()))
    overall = 100 * (grand[0] - grand[1]) / grand[0]
    with open(f"{out}/unify_flags.txt", "w") as fh:
        fh.write(f"grand_original={grand[0]}  grand_valid={grand[1]}  "
                 f"overall_doublet_pct={overall:.2f}\n")
        fh.write("\n".join(flags) if flags else "no anomalies\n")

    print(f"\n[summary] grand original={grand[0]:,}  valid={grand[1]:,}  "
          f"overall doublet rate={overall:.2f}%")
    print(f"[cluster17] original={bc.loc[17,'original_n']}  valid={bc.loc[17,'valid_n']}  "
          f"doublet%={bc.loc[17,'doublet_pct']}")
    if flags:
        print("[flags]")
        for f in flags:
            print("  -", f)
    print(f"\n[ok] wrote unified set to {out}/")


if __name__ == "__main__":
    main()
