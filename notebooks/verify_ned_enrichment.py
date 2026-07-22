# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "scipy"]
# ///
"""Reviewer 1 #3 (and #2) — the neuroendocrine states validate in independent human cohorts.

Reviewer-facing verification layer over the output of:
    scripts/integration_permanova_knn.py / annotate_validation_phenotypes.py
    -> data/integration_qc/validation_ne_enrichment.tsv  (also Supplementary Table T7)

This notebook RECOMPUTES the headline cross-cohort enrichment odds ratios live from the
deposited 2x2 counts, using the same Fisher exact test the pipeline uses, so the
Gao/Li odds ratios are derived in front of the reviewer rather than retyped.

Run interactively:        marimo edit notebooks/verify_ned_enrichment.py
Reproducible sandbox:     marimo edit --sandbox notebooks/verify_ned_enrichment.py
Export static HTML:       marimo export html notebooks/verify_ned_enrichment.py -o exports/verify_ned_enrichment.html
Run as a script:          python notebooks/verify_ned_enrichment.py
"""

import marimo

__generated_with = "0.23.4"
app = marimo.App(
    width="medium",
    app_title="Neuroendocrine states replicate in patient cohorts",
)


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _intro(mo):
    mo.md(r"""
    # Neuroendocrine states replicate in independent human cohorts

    *Addresses Reviewer 1, Comments 3 and 2.*

    > *"Do the neuroendocrine states identified in the LTL331/R model actually
    > correspond to neuroendocrine disease in independent patient data?"*

    ### Response in one sentence

    In both external cohorts, the cells our model assigns to the **neuroendocrine
    (NE)** program concentrate overwhelmingly in the epithelium that the cohort's
    own **clinical neuroendocrine-differentiation (NED)** call flags as
    neuroendocrine (Fisher odds ratio **134.6** in Gao and **31.4** in Li, both
    P far below 0.001). Because the cohorts do not fully co-embed (iLISI ≈ 1;
    see `verify_integration_qc`), this validation is at the label level, against
    the independent pathology call, not by a joint UMAP.

    Every odds ratio below is recomputed from the deposited 2x2 counts with the
    same Fisher exact test the pipeline uses.
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
    from scipy.stats import fisher_exact  # the pipeline's own enrichment test

    NED_TSV = repo_root / "data" / "integration_qc" / "validation_ne_enrichment.tsv"
    return NED_TSV, fisher_exact, np, pd


@app.cell(hide_code=True)
def _guard(NED_TSV, mo):
    # graceful message instead of a traceback if the input was not shipped / not yet built
    _missing = not NED_TSV.exists()
    if _missing:
        mo.callout(
            mo.md(
                f"""
                **Input not found:** `{NED_TSV}`

                Build it first (or fetch the deposited table):
                ```
                python scripts/integration_permanova_knn.py --adata <joint>.h5ad \\
                    --out data/integration_qc/
                ```
                """
            ),
            kind="warn",
        )
    mo.stop(_missing)
    return


@app.cell(hide_code=True)
def _recompute_header(mo):
    mo.md(r"""
    ## 1. Cross-cohort NE-vs-NED enrichment (recomputed)

    For each cohort we build the 2x2 table of *our* NE assignment against the
    cohort's *clinical* NED call and run a Fisher exact test. A high odds ratio
    means the two independent labels agree: the model's neuroendocrine cells are
    the clinically neuroendocrine cells.
    """)
    return


@app.cell
def _recompute(NED_TSV, fisher_exact, pd):
    dep = pd.read_csv(NED_TSV, sep="\t")
    rows = []
    for _, r in dep.iterrows():
        # 2x2: rows = our NE (yes/no), cols = clinical NED (yes/no)
        table = [[int(r["NE_NED"]), int(r["NE_nonNED"])],
                 [int(r["nonNE_NED"]), int(r["nonNE_nonNED"])]]
        odds, p = fisher_exact(table, alternative="two-sided")
        rows.append({
            "cohort": r["cohort"],
            "n_cells": sum(sum(x) for x in table),
            "OR_recomputed": round(float(odds), 2),
            "OR_deposited": round(float(r["odds_ratio"]), 2),
            "p_recomputed": p,
            "p_deposited": float(r["p_value"]),
            "matches": bool(abs(float(odds) - float(r["odds_ratio"])) < 0.5),
        })
    rec = pd.DataFrame(rows).set_index("cohort")
    all_match = bool(rec["matches"].all())
    or_gao = float(rec.loc["Gao", "OR_recomputed"]) if "Gao" in rec.index else float("nan")
    or_li = float(rec.loc["Li", "OR_recomputed"]) if "Li" in rec.index else float("nan")
    return all_match, or_gao, or_li, rec


@app.cell(hide_code=True)
def _verdict(all_match, mo, or_gao, or_li, rec):
    mo.callout(
        mo.md(
            f"""
            **Recomputed enrichment:** Gao odds ratio = **{or_gao:.1f}**,
            Li odds ratio = **{or_li:.1f}** — both far above 1 and matching the
            deposited values ({'all cohorts reproduce' if all_match else 'MISMATCH — check inputs'}).

            The model's neuroendocrine assignment tracks each cohort's independent
            clinical NED pathology call, so the neuroendocrine **endpoints** replicate
            in human disease even though the datasets do not fully co-embed.
            """
        ),
        kind="success" if all_match else "danger",
    )
    return


@app.cell
def _table(mo, rec):
    mo.ui.table(rec.reset_index(), page_size=10)
    return


@app.cell(hide_code=True)
def _provenance(mo):
    mo.md(r"""
    ---
    ### Provenance / full recompute

    The 2x2 counts read above come from the joint-integration validation:

    ```
    python scripts/integration_permanova_knn.py --adata <joint>.h5ad \
        --out data/integration_qc/        # -> validation_ne_enrichment.tsv
    ```

    This notebook re-derives the per-cohort Fisher odds ratios from those counts —
    lands on **main Fig. 6D + Supplementary Table T7 (R1 #3 / #2)**. The odds-ratio
    confidence intervals in the deposited table (Gao 60–300; Li 24–41) are reported
    in Fig. 6D.
    """)
    return


if __name__ == "__main__":
    app.run()
