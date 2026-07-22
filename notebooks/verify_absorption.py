# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "scipy", "matplotlib"]
# ///
"""Reviewer 1 #5 — cluster 17 is biased toward the ASCL1- NE fate.

Reviewer-facing verification layer over the output of:
    scripts/cluster17_absorption.py   (CellRank2 GPCCA absorption probabilities)

The heavy step (CellRank GPCCA on the velocity object, needs SLEPc) is NOT
re-run here. Instead this notebook RECOMPUTES the headline split from the
deposited per-cell fate-probability matrix, importing the *same* column-grouping
helper (`grp`) the pipeline uses — so the within-NE ASCL1+ fraction and the
Wilcoxon p-value are derived live, with only pandas/scipy.

Run interactively:        marimo edit notebooks/verify_absorption.py
Reproducible sandbox:     marimo edit --sandbox notebooks/verify_absorption.py
Export static HTML:        marimo export html notebooks/verify_absorption.py -o exports/verify_absorption.html
Run as a script:           python notebooks/verify_absorption.py
"""

import marimo

__generated_with = "0.23.4"
app = marimo.App(width="medium", app_title="Cluster 17 fate bias (ASCL1- vs ASCL1+)")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _intro(mo):
    mo.md(r"""
    # Cluster 17 fate bias (ASCL1- vs ASCL1+)

    *Addresses Reviewer 1, Comment 5.*

    > *"The claim that cluster 17 bifurcates into ASCL1+ and ASCL1- NEPC
    > states is asserted from marker expression. A quantitative,
    > statistically tested fate analysis (e.g. CellRank absorption
    > probabilities) is needed to support directional commitment."*

    ### Response in one sentence

    On this object's own CellRank terminal macrostates, cluster-17 cells
    carry a **median within-NE ASCL1+ fraction of 0.057** (i.e. ~0.94 toward
    ASCL1-), and a paired Wilcoxon test gives **P(ASCL1-) > P(ASCL1+),
    P = 4.76e-23** over the 128 cluster-17 cells.

    Recomputed below from the deposited fate-probability matrix by importing
    `grp` from `scripts/cluster17_absorption.py`.
    """)
    return


@app.cell(hide_code=True)
def _method_note(mo):
    mo.callout(
        mo.md(
            r"""
            **Which test is load-bearing.** The honest comparison is *within the
            NE compartment, forward terminals only*: P(ASCL1+) = terminal {10}
            vs P(ASCL1-) = terminals {1, 7, 9}. The NE-vs-non-NE contrast is
            **deliberately not** the headline — the non-NE "terminals" (0/4/6)
            are upstream clusters that GPCCA over-selects, so P(non-NE) is mostly
            backward-diffusion leakage from the symmetric connectivity kernel,
            not biology. We report the bifurcation test, not that one.
            """
        ),
        kind="warn",
    )
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
    from scipy.stats import wilcoxon
    from cluster17_absorption import grp  # the pipeline's own terminal-grouping helper

    ABS_DIR = repo_root / "data" / "forfig5" / "cluster17_absorption"
    BIF_DIR = repo_root / "data" / "forfig5" / "cr_bifurcation"
    # NE / non-NE / ASCL1+/- terminal definitions, identical to the script defaults
    NE = ["1", "7", "9", "10"]
    ASCL1_POS = ["10"]
    ASCL1_NEG = ["1", "7", "9"]
    return ABS_DIR, ASCL1_NEG, ASCL1_POS, BIF_DIR, grp, np, pd, wilcoxon


@app.cell(hide_code=True)
def _recompute_header(mo):
    mo.md(r"""
    ## 1. Recompute the cluster-17 fate split

    Per-cell fate probabilities are loaded from
    `fate_probabilities_allcells.tsv`; the barcode→cluster map (to isolate
    the 128 cluster-17 cells) comes from `fate_probabilities_seedmain.tsv`,
    which carries the aligned `seurat_clusters` column.
    """)
    return


@app.cell
def _recompute(ABS_DIR, ASCL1_NEG, ASCL1_POS, BIF_DIR, grp, np, pd, wilcoxon):
    fp = pd.read_csv(ABS_DIR / "fate_probabilities_allcells.tsv", sep="\t", index_col=0)
    fp.columns = [str(c) for c in fp.columns]
    seed = pd.read_csv(BIF_DIR / "fate_probabilities_seedmain.tsv", sep="\t", index_col=0)
    cl = seed["seurat_clusters"].astype(str)

    cl17_bc = cl.index[cl == "17"]
    sub = fp.loc[fp.index.intersection(cl17_bc)]
    n_src = len(sub)

    p_pos, _ = grp(sub, ASCL1_POS)
    p_neg, _ = grp(sub, ASCL1_NEG)
    ne_tot = (p_pos + p_neg).replace(0, np.nan)
    frac_pos = p_pos / ne_tot

    median_frac_pos = float(np.nanmedian(frac_pos))
    stat, pval = wilcoxon(p_neg.values, p_pos.values, alternative="greater")
    return frac_pos, median_frac_pos, n_src, p_neg, p_pos, pval


@app.cell(hide_code=True)
def _recompute_verdict(median_frac_pos, mo, n_src, pval):
    mo.callout(
        mo.md(
            f"""
            **Cluster 17 (n = {n_src}) — recomputed:**
            median within-NE **P(ASCL1+) = {median_frac_pos:.4f}**
            ⇒ ~**{1 - median_frac_pos:.2f}** toward **ASCL1-**.

            Paired Wilcoxon, P(ASCL1-) > P(ASCL1+): **P = {pval:.2e}**.

            Cluster 17 commits predominantly to the **ASCL1- NE fate**.
            """
        ),
        kind="success",
    )
    return


@app.cell
def _fate_plot(frac_pos, np, p_neg, p_pos):
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2))

    # left: per-cell ASCL1- vs ASCL1+ absorption mass
    ax1.scatter(p_pos.values, p_neg.values, s=14, c="#7570b3", edgecolor="k", lw=0.2)
    lim = [0, max(p_pos.max(), p_neg.max()) * 1.05]
    ax1.plot(lim, lim, ls=":", c="grey", lw=0.8)
    ax1.set_xlim(lim); ax1.set_ylim(lim)
    ax1.set_xlabel("P(ASCL1+)  terminal {10}")
    ax1.set_ylabel("P(ASCL1-)  terminals {1,7,9}")
    ax1.set_title("cluster-17 cells: ASCL1- dominates")

    # right: within-NE ASCL1+ fraction distribution
    ax2.hist(frac_pos.dropna().values, bins=20, color="#d95f02", edgecolor="k", lw=0.3)
    ax2.axvline(np.nanmedian(frac_pos), ls="--", c="k", lw=1,
                label=f"median {np.nanmedian(frac_pos):.3f}")
    ax2.set_xlabel("within-NE fraction ASCL1+")
    ax2.set_ylabel("cluster-17 cells")
    ax2.set_title("mass piled near 0 (ASCL1-)")
    ax2.legend(frameon=False)
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _deposited_header(mo):
    mo.md(r"""
    ## 2. Deposited figure & summary table
    """)
    return


@app.cell
def _deposited_image(BIF_DIR, mo):
    img = BIF_DIR / "intermediate_fate_violin.png"
    out = (mo.image(str(img), alt="cluster-17 fate violin")
           if img.exists() else mo.md(f"_(missing image: `{img}`)_"))
    out
    return


@app.cell
def _summary_table(ABS_DIR, mo, pd):
    summ = pd.read_csv(ABS_DIR / "cl17_fate_summary.tsv", sep="\t")
    mo.ui.table(summ.round(6), page_size=5)
    return


@app.cell(hide_code=True)
def _provenance(mo):
    mo.md(r"""
    ---
    ### Provenance / full recompute

    The fate-probability matrix is produced by the CellRank2 GPCCA run
    (needs the velocity object + SLEPc — not re-run in this notebook):

    ```
    python scripts/cluster17_absorption.py \
        --h5ad ltl331_velocity_cellrank.h5ad \
        --source 17 --term-key term_states_fwd \
        --ne 1 7 9 10 --non-ne 0 4_1 4_2 6 \
        --ascl1-pos 10 --ascl1-neg 1 7 9 \
        --outdir data/forfig5/cluster17_absorption
    ```

    This notebook re-derives the within-NE ASCL1± split and the Wilcoxon
    p-value from the deposited per-cell matrix using the pipeline's own
    `grp` helper — lands on **Fig 4F + Supp S14 (R1 #5, R2 #6)**.
    """)
    return


if __name__ == "__main__":
    app.run()
