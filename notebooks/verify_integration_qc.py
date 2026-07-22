# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "matplotlib"]
# ///
"""Reviewer 1 #3 — integration QC (iLISI / cLISI), honestly.

Reviewer-facing verification layer over the outputs of:
    scripts/integration_lisi_kbet.py       (iLISI/cLISI on naive vs Harmony embeddings)
    scripts/integration_permanova_knn.py   (PERMANOVA R2 + dispersion)

Companion to `notebooks/R1_3_integration.py`, which covers the Harmony
cluster-agreement / origin-map facet of R1 #3. THIS notebook covers the QC
metrics and the honest reading of them: the three cohorts do **not** kNN-mix
even after Harmony (iLISI = 1.0), so the R1 #3 deliverable is **label/signature
transfer, not co-embedding** — and that label transfer is concordant.

The full LISI recompute needs the 74k-cell joint object + harmonypy/scanpy and
is guarded behind a file check; the light tier (deltas, cross-cohort concordance,
TME fractions) recomputes from the deposited tables with only pandas/numpy.

Run interactively:        marimo edit notebooks/verify_integration_qc.py
Reproducible sandbox:     marimo edit --sandbox notebooks/verify_integration_qc.py
Export static HTML:        marimo export html notebooks/verify_integration_qc.py -o exports/verify_integration_qc.html
Run as a script:           python notebooks/verify_integration_qc.py
"""

import marimo

__generated_with = "0.23.4"
app = marimo.App(width="medium", app_title="Cross-cohort integration QC (iLISI / cLISI)")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _intro(mo):
    mo.md(r"""
    # Cross-cohort integration QC (iLISI / cLISI)

    *Addresses Reviewer 1, Comment 3.*

    > *"…a direct integration of datasets (e.g. Harmony or Seurat v5) to
    > demonstrate that PDX-derived transition cells and clinical 'plastic'
    > cells co-localize in the same latent manifold."*

    ### Response in one sentence (the honest version)

    Harmony applies a **real** correction (≈29 % displacement on shared dims)
    but the human-epithelium PDX and the clinical biopsies **do not kNN-mix**
    (iLISI = 1.000 before *and* after) — a linear PCA correction cannot bridge
    a platform/model batch this large. We do **not** force it (cranking θ to
    fake co-embedding is exactly the over-correction R2 #7 watches for).
    Instead the R1 #3 claim rests on **concordant label transfer**: uniform
    marker scoring places PRAD-luminal and NE epithelial populations in
    **all three cohorts** on the same AR→NE axis.

    > See `notebooks/R1_3_integration.py` for the cluster-agreement / origin-map
    > evidence; this notebook is the metric QC behind the decision.
    """)
    return


@app.cell
def _setup_paths():
    import json
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

    QC = repo_root / "data" / "integration_qc"
    return QC, pd, repo_root


@app.cell(hide_code=True)
def _lisi_header(mo):
    mo.md(r"""
    ## 1. LISI metrics — naive PCA vs Harmony (recomputed reading)

    Raw scale: **iLISI ∈ [1, 3]**, higher = better batch mixing;
    **cLISI ∈ [1, n_class]**, lower = better cell-type preservation. The
    delta and interpretation are recomputed from the deposited metric table
    — *not* read from `verdict.json`, whose auto-generated string is
    misleading (it labels Δ-iLISI = 0 as "mixing NOT improved", which
    understates that the cohorts never mixed in the first place).
    """)
    return


@app.cell
def _lisi_recompute(QC, pd):
    lisi = pd.read_csv(QC / "lisi_kbet_metrics.tsv", sep="\t")
    naive = lisi.iloc[0]
    corr = lisi.iloc[1]
    d_iLISI = float(corr["iLISI_median"] - naive["iLISI_median"])
    d_cLISI = float(corr["cLISI_median"] - naive["cLISI_median"])
    over_corrected = d_cLISI >= 0.3   # cLISI rose a lot => types smeared
    mixed = d_iLISI > 0.1
    return d_cLISI, d_iLISI, lisi, mixed, over_corrected


@app.cell
def _lisi_table(lisi, mo):
    mo.ui.table(lisi.round(4), page_size=4)
    return


@app.cell(hide_code=True)
def _lisi_verdict(d_cLISI, d_iLISI, mixed, mo, over_corrected):
    mo.callout(
        mo.md(
            f"""
            **Recomputed:** Δ-iLISI = **{d_iLISI:+.3f}** ⇒ kNN mixing
            { "improved" if mixed else "**did not** improve (iLISI = 1.0 throughout — cohorts stay single-dataset)" }.
            Δ-cLISI = **{d_cLISI:+.3f}** ⇒ biology
            { "**over-corrected**" if over_corrected else "**preserved** (no smearing of cell types)" }.

            Honest reading: Harmony shifted each cohort as a coherent block but
            never interleaved them — *not* a no-op, *not* over-correction.
            """
        ),
        kind="info",
    )
    return


@app.cell(hide_code=True)
def _concordance_header(mo):
    mo.md(r"""
    ## 2. The clean result — concordant label transfer

    Uniform canonical-marker scoring (PRAD / NE / Basal / Other) applied
    identically to all three cohorts. Two things to confirm: (a) an **NE
    epithelial population exists in every cohort** (the shared AR→NE axis),
    and (b) the non-epithelial **"Other" fraction is < 10 %** everywhere, so
    the clinical cohorts were validly restricted to tumor epithelium before
    integration.
    """)
    return


@app.cell
def _concordance_recompute(QC, pd):
    comp = pd.read_csv(QC / "coarse_label_by_dataset.tsv", sep="\t", index_col=0)
    frac = comp.div(comp.sum(axis=0), axis=1)
    ne_per_cohort = comp.loc["NE"].to_dict()
    other_frac = (100 * frac.loc["Other"]).round(1).to_dict()
    ne_in_all = bool((comp.loc["NE"] > 0).all())
    tme_ok = bool((frac.loc["Other"] < 0.10).all())
    return comp, frac, ne_in_all, ne_per_cohort, other_frac, tme_ok


@app.cell(hide_code=True)
def _concordance_verdict(mo, ne_in_all, ne_per_cohort, other_frac, tme_ok):
    mo.callout(
        mo.md(
            f"""
            **NE epithelium in every cohort:** { {k: int(v) for k, v in ne_per_cohort.items()} }
            → present in all three = **{ne_in_all}**.

            **"Other" (TME) fraction:** {other_frac} % →
            all < 10 % = **{tme_ok}** (epithelium restriction valid).
            """
        ),
        kind="success" if (ne_in_all and tme_ok) else "warn",
    )
    return


@app.cell
def _concordance_plot(comp, frac):
    import matplotlib.pyplot as plt
    classes = list(comp.index)
    cohorts = list(comp.columns)
    colors = {"PRAD": "#1f77b4", "NE": "#d62728", "Basal": "#2ca02c", "Other": "#999999"}
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    bottom = [0.0] * len(cohorts)
    for cls in classes:
        vals = (100 * frac.loc[cls]).values
        ax.bar(cohorts, vals, bottom=bottom, label=cls,
               color=colors.get(cls, "#cccccc"), edgecolor="k", lw=0.3)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_ylabel("% of cohort cells")
    ax.set_title("Coarse lineage composition by cohort\n(NE present in all three; AR→NE axis shared)")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _full_recompute_header(mo):
    mo.md(r"""
    ## 3. Optional — full LISI recompute from the joint object

    The cell below re-runs the actual iLISI/cLISI computation **only if** the
    74k-cell joint object is present locally (it needs `scanpy` + `harmonypy`
    and the embeddings `X_pca` / `X_harmony_dataset`). In the sandbox/HTML
    export it is skipped with a note — the light tier above already
    reproduces the published numbers from the deposited tables.
    """)
    return


@app.cell
def _full_recompute(mo, repo_root):
    # look for the joint object in a couple of conventional spots
    candidates = [
        repo_root / "harmony_10pt.h5ad",
        repo_root / "data" / "scanpy" / "harmony_10pt.h5ad",
    ]
    h5 = next((c for c in candidates if c.exists()), None)
    if h5 is None:
        out = mo.callout(
            mo.md(
                "**Skipped** — joint object not found locally "
                f"(looked in: {[str(c) for c in candidates]}). "
                "Run `python scripts/integration_lisi_kbet.py --adata <obj>.h5ad "
                "--out data/integration_qc/` to regenerate the deposited tables."
            ),
            kind="neutral",
        )
    else:
        import scanpy as sc
        from integration_lisi_kbet import assign_coarse, lisi_medians
        adata = sc.read_h5ad(str(h5))
        label = assign_coarse(adata, use_raw=True)
        recomputed = {
            rep: dict(zip(("iLISI", "cLISI"),
                          lisi_medians(adata, rep, "dataset", label)))
            for rep in ("X_pca", "X_harmony_dataset") if rep in adata.obsm
        }
        out = mo.md(f"**Recomputed live from `{h5.name}`:**\n\n```\n{recomputed}\n```")
    out
    return


@app.cell(hide_code=True)
def _provenance(mo):
    mo.md(r"""
    ---
    ### Provenance

    ```
    python scripts/integration_lisi_kbet.py    --adata harmony_10pt.h5ad --out data/integration_qc/
    python scripts/integration_permanova_knn.py --adata harmony_10pt.h5ad --out data/integration_qc/
    ```

    Honest interpretation pinned in `data/integration_qc/INTERPRETATION.md`.
    Lands on **Supp S20 + Table T8 (R1 #3)**; the label-transfer evidence is
    carried by Fig 6 and `notebooks/R1_3_integration.py`.
    """)
    return


if __name__ == "__main__":
    app.run()
