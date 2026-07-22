# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "scipy", "matplotlib"]
# ///
"""Reviewer 1 #1 — cluster 17 is not a doublet artefact.

Reviewer-facing verification layer over the outputs of:
    scripts/doublet_unify.py        (DoubletFinder original vs valid -> unified tables)
    scripts/doublet_robustness.py   (per-cluster Fisher enrichment + composition stats)

This notebook RECOMPUTES the headline numbers live from the deposited unified
tables, importing the *same* helper (`bh`) the pipeline uses — so the cluster-17
odds ratio and enrichment call are derived in front of the reviewer, not retyped.

Run interactively:        marimo edit notebooks/verify_doublets.py
Reproducible sandbox:     marimo edit --sandbox notebooks/verify_doublets.py
Export static HTML:        marimo export html notebooks/verify_doublets.py -o exports/verify_doublets.html
Run as a script:           python notebooks/verify_doublets.py
"""

import marimo

__generated_with = "0.23.4"
app = marimo.App(
    width="medium",
    app_title="Cluster 17 is not a doublet artefact",
)


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _intro(mo):
    mo.md(r"""
    # Cluster 17 is not a doublet artefact

    *Addresses Reviewer 1, Comment 1.*

    > *"The authors should rule out that the rare intermediate population
    > (cluster 17) is a technical doublet artefact rather than a genuine
    > transcriptional state."*

    ### Response in one sentence

    DoubletFinder flags cluster 17 as **depleted** of doublets (rate 1.56 %
    vs 4.17 % overall; Fisher OR = **0.36**, not significant), and removing
    every called doublet leaves cluster composition (r = 0.999) and cluster
    17's wk12–16 temporal peak essentially unchanged.

    Every number below is recomputed from the deposited tables by importing
    `bh` (Benjamini–Hochberg) from `scripts/doublet_robustness.py`.
    """)
    return


@app.cell
def _setup_paths():
    import sys
    from pathlib import Path
    try:
        repo_root = Path(__file__).resolve().parent.parent
    except NameError:  # interactive
        repo_root = Path.cwd().resolve()
    _paths = [str(repo_root), str(repo_root / "scripts")]
    sys.path[:0] = [q for q in _paths if q not in sys.path]

    import numpy as np
    import pandas as pd
    from scipy.stats import fisher_exact
    from doublet_robustness import bh  # the pipeline's own BH helper

    UNIFIED = repo_root / "data" / "doubletfinder" / "unified"
    return UNIFIED, bh, fisher_exact, np, pd


@app.cell(hide_code=True)
def _enrich_header(mo):
    mo.md(r"""
    ## 1. Per-cluster doublet enrichment (recomputed)

    Fisher exact test, doublets-in-cluster vs rest, BH-corrected. Cluster 17
    should appear with **OR < 1** and **not** in the enriched set (the
    enriched clusters are the cycling "engine" states 3/5/6/8/9/10/13).
    """)
    return


@app.cell
def _enrich_recompute(UNIFIED, bh, fisher_exact, pd):
    bc = (pd.read_csv(UNIFIED / "doublet_by_cluster.tsv", sep="\t")
          .set_index("cluster").sort_index())
    tot_o, tot_v = int(bc["original_n"].sum()), int(bc["valid_n"].sum())
    tot_d = tot_o - tot_v
    overall = 100 * tot_d / tot_o

    rows = []
    for c in bc.index:
        d_c = max(int(bc.loc[c, "doublets"]), 0)
        keep_c = int(bc.loc[c, "valid_n"])
        odds, p = fisher_exact([[d_c, keep_c], [tot_d - d_c, tot_v - keep_c]])
        rows.append({"cluster": c, "doublet_pct": float(bc.loc[c, "doublet_pct"]),
                     "odds_ratio": round(float(odds), 3), "fisher_p": p})
    enr = pd.DataFrame(rows).set_index("cluster")
    enr["bh_q"] = bh(enr["fisher_p"].values)
    enr["enriched"] = (enr["odds_ratio"] > 1) & (enr["bh_q"] < 0.05)

    cl17_or = float(enr.loc[17, "odds_ratio"])
    cl17_enriched = bool(enr.loc[17, "enriched"])
    enriched_clusters = [int(c) for c in enr.index[enr["enriched"]]]
    return (
        bc,
        cl17_enriched,
        cl17_or,
        enr,
        enriched_clusters,
        overall,
        tot_o,
        tot_v,
    )


@app.cell(hide_code=True)
def _enrich_verdict(cl17_enriched, cl17_or, enriched_clusters, mo, overall):
    mo.callout(
        mo.md(
            f"""
            **Cluster 17 — recomputed:** doublet OR = **{cl17_or:.3f}**,
            enriched = **{cl17_enriched}** (overall doublet rate {overall:.2f} %).

            **Enriched clusters (BH < 0.05):** {enriched_clusters} — the cycling
            engines, *not* cluster 17. Cluster 17's intermediate state is
            doublet-**depleted**, the opposite of a doublet artefact.
            """
        ),
        kind="success" if not cl17_enriched else "danger",
    )
    return


@app.cell
def _enrich_table(enr, mo):
    mo.ui.table(enr.reset_index().round(4), page_size=18)
    return


@app.cell(hide_code=True)
def _comp_header(mo):
    mo.md(r"""
    ## 2. Composition & trajectory preserved after doublet removal

    If cluster 17 were doublet-driven, deleting the called doublets would
    distort the cluster fractions and collapse cluster 17's temporal peak.
    It does neither.
    """)
    return


@app.cell
def _comp_recompute(bc, np, tot_o, tot_v):
    o_frac = bc["original_n"] / tot_o
    v_frac = bc["valid_n"] / tot_v
    pearson = float(np.corrcoef(bc["original_n"], bc["valid_n"])[0, 1])
    cosine = float(np.dot(o_frac, v_frac) /
                   (np.linalg.norm(o_frac) * np.linalg.norm(v_frac)))
    max_shift_pp = float((100 * (v_frac - o_frac)).abs().max())
    max_shift_cluster = int((v_frac - o_frac).abs().idxmax())
    return cosine, max_shift_cluster, max_shift_pp, pearson


@app.cell
def _comp_temporal(UNIFIED, np, pd):
    long = pd.read_csv(UNIFIED / "doublet_counts_long.tsv", sep="\t")
    piv_o = long.pivot_table(index="timepoint", columns="cluster",
                             values="original_n", aggfunc="sum").fillna(0)
    piv_v = long.pivot_table(index="timepoint", columns="cluster",
                             values="valid_n", aggfunc="sum").fillna(0)
    tps = long.sort_values("tp_order")["timepoint"].drop_duplicates().tolist()
    piv_o, piv_v = piv_o.reindex(tps), piv_v.reindex(tps)
    of = piv_o.div(piv_o.sum(1), axis=0)
    vf = piv_v.div(piv_v.sum(1), axis=0)
    per_tp_r = {t: float(np.corrcoef(of.loc[t], vf.loc[t])[0, 1]) for t in tps}
    cl17_traj = pd.DataFrame({"original": 100 * of[17], "valid": 100 * vf[17]})
    min_tp_r = min(per_tp_r.values())
    return cl17_traj, min_tp_r, tps


@app.cell(hide_code=True)
def _comp_verdict(
    cosine,
    max_shift_cluster,
    max_shift_pp,
    min_tp_r,
    mo,
    pearson,
):
    mo.callout(
        mo.md(
            f"""
            **Composition preserved:** Pearson r = **{pearson:.4f}**, cosine =
            **{cosine:.5f}**; largest per-cluster fraction shift = {max_shift_pp:.2f}
            pp (at cluster {max_shift_cluster}).
            **Trajectory preserved:** worst per-timepoint composition correlation
            = **{min_tp_r:.4f}**.
            """
        ),
        kind="info",
    )
    return


@app.cell
def _cl17_traj_plot(cl17_traj, tps):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(tps, cl17_traj["original"], "-o", ms=4, c="#3182bd", label="original")
    ax.plot(tps, cl17_traj["valid"], "--s", ms=4, c="#e7298a", label="doublets removed")
    ax.set_xlabel("timepoint")
    ax.set_ylabel("cluster 17 (% of timepoint)")
    ax.set_title("Cluster 17 wk12–16 peak survives doublet removal")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _figure_header(mo):
    mo.md(r"""
    ## 3. Full supplementary panel (deposited)

    The 5-panel figure emitted by `doublet_robustness.py` (the version that
    ships as a supplementary figure).
    """)
    return


@app.cell
def _figure_image(UNIFIED, mo):
    img = UNIFIED / "doublet_robustness.png"
    out = (mo.image(str(img), alt="Doublet robustness 5-panel")
           if img.exists() else mo.md(f"_(missing image: `{img}`)_"))
    out
    return


@app.cell(hide_code=True)
def _provenance(mo):
    mo.md(r"""
    ---
    ### Provenance / full recompute

    The unified tables read above come from the upstream DoubletFinder run:

    ```
    python scripts/doublet_unify.py        # original vs valid -> data/doubletfinder/unified/
    python scripts/doublet_robustness.py \
        --unified data/doubletfinder/unified --target-cluster 17
    ```

    This notebook re-derives the cluster-17 OR, enrichment call, composition
    and temporal-preservation statistics from those tables using the
    pipeline's own `bh` helper — lands on **Supp S6 + Table T5 (R1 #1)**.
    """)
    return


if __name__ == "__main__":
    app.run()
