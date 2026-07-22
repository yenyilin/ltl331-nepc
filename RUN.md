# RUN ORDER — reproducible workflow (LTL331/R NEPC)

Canonical execution order. **Stages reflect data dependencies, not figure numbers**
(scripts stay flat under `scripts/python/`; this file is the run order,
`FIGURES.md` is the figure→script lookup). A later stage consumes objects/labels built
by earlier ones.

Note on "integration": there are two. The **within-LTL331 assembly** (timepoints + ERCC
batch) is preprocessing — stage 00. The **cross-cohort Harmony** (LTL331 + Gao + Li) is the
*clinical-validation* integration and is **late** — stage 05 — because it transfers the
established LTL331 trajectory/signatures onto clinical cells and needs tumor-epithelium
selection (inferCNV) first.

Envs: `cci` (anndata 0.12.10) and `cellrank` (anndata 0.10.7). h5ad cross-env rule: write in
`cci`, read in `cellrank`; `make_portable_h5ad.py` bridges (velocity object: use the original,
which already reads in 0.10.7). See `env/` for exact versions.

---

## 00 · Preprocess & QC  → the LTL331 object
- per-sample CellRanger aggr → QC; **within-LTL331 assembly** (8 timepoints, ERCC batch boundary)
- `scripts/python/doublet_unify.py` → `doublet_robustness.py`   (R1#1 doublet sensitivity)
  - in: `data/doubletfinder/counts_{original,valid}.tsv` (before/after count matrices; see that dir's README)
  - out: `unified/` tables + `doublet_robustness.{pdf,png}` = **Fig. S1** (does not touch the h5ad)
- `scripts/python/canonicalize_cellid.py`  (stable join key across objects)
- out: LTL331 base h5ad (lognorm in .raw)

## 01 · Clustering & cluster taxonomy  (LTL331 core)
- Seurat clustering → 18-cluster taxonomy, `seurat_clusters`, markers (DEG.txt)
- `scripts/python/make_signatures.py`  (DEG → per-cluster signatures; feeds 03 & 05)
- canonical `data/cluster_order.json` (AR→NE phenotypic axis), `cluster_colors_18.json`
- out: clustered LTL331 object + signatures + cluster order

## 02 · Trajectory  (parallel branch off 01)
- `scripts/python/velocity_dynamical.py` (scVelo) → `cellrank_bifurcation.py` (R1#5 ASCL1±)
- transition flux / pseudotime QC
- out: velocity/CellRank object (use the ORIGINAL h5ad in cellrank env; reads in 0.10.7)

## 03 · Regulons  (parallel branch off 01)
- `scripts/python/pyscenic_run.py` (SCENIC + NKX2-1/ONECUT2/BRN2; R1#4)
- `scripts/python/decoupler_tf_activity.py --pseudobulk-by seurat_clusters` (LTL331)
- `scripts/python/scenic_decoupler_consensus.py`  (method-orthogonal regulon confirmation)
- out: regulon activity + SCENIC↔decoupleR consensus

## 04 · CNV  (parallel branch off 01; prerequisite for 05 clinical tumor selection)
- inferCNV per timepoint/cluster (computed in R, `infercnv`; versions in `config/params.yaml`) -> CNV calls deposited and consumed by the Python plotting scripts (R2#2)
- Gao/Li tumor-epithelium selection (inferCNV + markers) — gates entry to stage 05
- out: PDX tumor/clonal confirmation; clinical tumor-cell whitelist

## 05 · Clinical validation  (cross-cohort — LATE; depends on 01,03,04)
- cross-cohort **Harmony** integration (LTL331 + Gao + Li, tumor epithelium)
- `scripts/python/integration_lisi_kbet.py` + `integration_permanova_knn.py`  (R1#3 integration QC)
- `scripts/python/harmony_cluster_agreement.py --use-existing-cols`  (resolution: Louvain 0.55)
- `scripts/python/compare_resolution_trajectory.py`, `plot_origin_vs_harmony.py`
  (`--cluster-order` + `--phenotype-col`), `cluster_origin_correspondence.py`
- `scripts/python/annotate_validation_phenotypes.py --report`  (NED→NE enrichment)
- `scripts/python/clinical_signature_validation.py`  (patient × disease-axis signature transfer)
- `scripts/python/decoupler_tf_activity.py --pseudobulk-by dataset,ned_status`  (cross-cohort regulons)
- out: Fig 6 integration map + clinical-validation supplements

## 06 · Figures
- per-figure assembly from all stages; see `FIGURES.md` for panel→script. Uses
  `cluster_order.json` / `cluster_colors_18.json` for consistent ordering/palette.

## 07 · Reviewer-facing verification notebooks  (marimo; recompute, don't reprint)
Each notebook re-derives a Tier-1 headline live from the deposited tables by importing
the pipeline's own helper, so the number is computed in front of the reviewer. PEP 723
headers → `uv` provisions an isolated env, sidestepping the main analysis env above.
- `notebooks/verify_doublets.py`        (R1#1) — imports `bh`; recomputes cl17 OR=0.364, not enriched; comp r=0.999
- `notebooks/verify_absorption.py`      (R1#5) — imports `grp`; recomputes within-NE ASCL1+ 0.057, P=4.76e-23 (n=128)
- `notebooks/verify_integration_qc.py`  (R1#3) — recomputes iLISI Δ=0 / cLISI Δ=0.02 + cross-cohort label concordance
- `notebooks/verify_ned_enrichment.py`  (R1#3/#2) — recomputes NE-vs-clinical-NED Fisher OR (Gao 134.6 / Li 31.4) from the deposited 2×2 (= Table T7)

```
marimo edit --sandbox notebooks/verify_doublets.py            # interactive, isolated env
marimo export html --sandbox notebooks/verify_doublets.py \
    -o exports/verify_doublets.html                           # static HTML for the response letter
```

---

### Stage → figure (capstone view)
- Figs 2–3  ← stages 00–01 (LTL331 taxonomy/composition)
- Fig 4     ← stage 02 (trajectory)
- Fig 5     ← stage 03 (regulons)
- **Fig 6  ← stage 05 (clinical validation)** — late, consistent with integration being late
- Supp CNV ← stage 04
