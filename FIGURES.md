# Figure → Script Map

Every main and supplementary figure panel and the script that produces it.

Legend: **[HAVE]** = dedicated script in this repository; **[orig]** = produced by the
original analysis pipeline; **[exp]** = commercial-software output, no reproducible script.
---

## Main figures

| Fig | Panel | Content | Script | Status |
|-----|-------|---------|--------|--------|
| 1 | A | PSA + tumor volume time course | - | orig (measurement table, not scRNA) |
| 1 | B | H&E + IHC (AR/PSA/NCAM1) | — | exp (no code) |
| 2 | A | UMAP of 51,726 cells, 18 clusters | `scripts/python/plot_umap_panels.py --color-by seurat_clusters` | orig / ~ |
| 2 | B | **same UMAP colored by timepoint** | `scripts/python/plot_timepoint_panel.py` | [HAVE] |
| 2 | C | AR & NCAM1 expression on UMAP | `scripts/python/plot_umap_panels.py --which continuous` (or orig feature plots) | ~ |
| 2 | D | dot plot, canonical lineage markers per cluster (**+ luminal/basal cytokeratins**), adeno→NE ordered | `scripts/python/plot_marker_dotplot.py` (+ `data/fig2c_main.tsv`) | [HAVE] |
| 2 | E | per-cluster PRAD/AR-pathway (blue) & neuroendocrine (red) **enrichment scores** | `scripts/python/plot_enrichment_violin.py` | [HAVE] |
| 2 | F | **basal × EMT score per cell** (cl17 highlighted; basal-*like*, not classical p63⁺/KRT5⁺) | `scripts/python/plot_cl17_basal_vs_emt.py` (+ `data/cl17_basal_emt_signatures.tsv`) | [HAVE] |
| 3 | A | bar: cluster proportions across timepoints | `scripts/python/plot_cluster_proportions.py` | [HAVE] |
| 3 | B | ridge: cells/timepoint along the AR→NE axis with AR-high/AR-low/Intermediate/NEPC stage bands | `scripts/python/plot_fig3b_timepoint_ridge.py` | [HAVE] |
| 3 | C | heatmap Hallmark GSEA NES *(was 3D)* | - | [orig] |
| 3 | D | **cluster 4/12 temporal dynamics** — within-timepoint % over 8 tps, cl4 survival-peak 58%@wk16 then lost as NE→99%; cl12 transient | `scripts/python/plot_fig3e_cl4_12_dynamics.py` → `figures/fig3/fig3E_cl4_12_dynamics.{pdf,png}` | [DONE] |
| 3 | E | cluster-17 EMT/Progenitor/Stemness pathways (C2:CGP) | - | [orig] |
| 3 | F | violin proliferation scores | `scripts/python/plot_fig3f_proliferation_violin.py` | [HAVE] |
| 3 | — | KRAS up/down disambiguation → Table T6 + text | `scripts/python/pathway_revision.py` | [HAVE] (no GSEA re-run) |
| 4 | A | PAGA connectivity (edge thr 0.09; `--flip` for orientation) | `scripts/python/paga_trajectory.py`; **`verify_paga_obligate.py`** (cut-vertex / obligate-gateway test: max PRAD↔NE conn 0.036, cl17↔cl13 1.0 — Results + Methods, R2-minor-#3) | [HAVE] |
| 4 | B | RNA-velocity stream embedding (dynamical) | `scripts/python/velocity_dynamical.py` | [HAVE] |
| 4 | C | velocity pseudotime | `scripts/python/velocity_dynamical.py` | [HAVE] |
| 4 | D | **two NE terminal states** (cl7 ASCL1−, cl10 ASCL1+) on UMAP | compute `scripts/python/cellrank_bifurcation.py` ⚙ → render `scripts/python/plot_fate_umap.py --terminals-panel` | [HAVE] |
| 4 | E | CellRank 2 fate probabilities on UMAP (within-NE ASCL1−/+ ratio) | `scripts/python/plot_fate_umap.py --within-ne-ratio` | [HAVE] |
| 4 | F | cluster-17 absorption-prob violin + test (ASCL1− 0.94 / P=4.76e-23) | render `scripts/python/plot_fate_umap.py` · compute `scripts/python/cluster17_absorption.py` + `scripts/python/cluster17_transition_flux.py` ⚙ | [HAVE] |
| 5 | A | SCENIC regulon map (size=RSS, color=AUCell; canonical NKX2-1/ONECUT2/POU3F2 grayed "not recovered by SCENIC") | `scripts/python/plot_scenic_regulon_map.py` (RSS) + `plot_scenic_auc_map.py` (AUCell) | [HAVE] |
| 5 | B | decoupleR TF activity (ASCL1 concordant with SCENIC at cl10; POU3F2/INSM1 active at cl10; MSX1 positive at the gateway; **`--zscore-pseudobulk`**; reindex cols→cluster_order) + ONECUT2 expr inset | `scripts/python/plot_tf_activity_heatmap.py`; `tf_expression_by_cluster.py` (ONECUT2) | [HAVE] |
| 5 | C | TF dynamics along pseudotime, both branches (17→13→10 / →7) | `scripts/python/plot_branch_tf_tracks.py --branches data/branches_conservative.json` | [HAVE] |
| 5 | D | **chromatin anchor** — clinical NEPC ATAC footprints (Shrestha; ASCL1/NEUROD1/POU3F2 in AR−/NE+); legend: clinical corroboration, not matched ATAC. **Run `--curated`** (NE + EMT/bridge story rows); stem `fig5D`. Full 17-TF atlas (no flag) → S15A | `scripts/python/plot_shrestha_chromatin.py --curated` | [HAVE] |
| 6 | A | joint Harmony UMAP, faceted by cohort × coarse lineage (orientation only) | `scripts/python/plot_umap_panels.py` | [HAVE] |
| 6 | B | MetaNeighbor cross-cohort AUROC (NE ASCL1+ 0.83, adeno 0.82; ASCL1− 0.63, basal/other ≈chance = internal controls) | `scripts/python/plot_terminal_metaneighbor.py` | [HAVE] |
| 6 | C | reference-projection label transfer (adeno→adeno 94%, NE→NE 57%; NE lands ASCL1−) | `scripts/python/plot_terminal_projection.py` | [HAVE] |
| 6 | D | NED→NE enrichment forest (Fisher OR 134.6 Gao / 31.4 Li, P<1e-10) | `scripts/python/ned_ne_enrichment.py` | [HAVE] |


---

## Supplementary figures (S1–S26)

Pipeline is **Python** (`scripts/python/*.py`). Status: **[HAVE]** dedicated script ·
**~** shared script · **[orig]** original-pipeline analysis · **[exp]** commercial software output.

| Supp | Content | Script | Status |
|------|---------|--------|--------|
| S1 | **Doublet sensitivity** (cl17 OR≈0.36; r=0.999) | `doublet_unify.py` → `doublet_robustness.py` | [HAVE] |
| S2 | Parental (927) + LTL331 CN profile (Nexus) | `nexus_cnv_loader.py` + `plot_nexus_cnv_heatmap.py` / `plot_nexus_cnv_curves.py` | [exp]/[HAVE] |
| S3 | WES CN profiles (preCx/16/20/22) | `nexus_cnv_loader.py` (Nexus output) | [orig] |
| S4 | Sample-level inferCNV (by timepoint + ungrouped) | inferCNV run (R pkg) → `infercnv_loader.py` | ~ |
| S5 | **Cluster-level inferCNV — cl17 malignant** | `infercnv_loader.py` + `infercnv_clonal_check.py` / `wk16_cnv_per_cluster.py` | ~ |
| S6 | **WES↔inferCNV per-chromosome concordance** | `compare_infercnv_vs_wes.py` | [HAVE] |
| S7 | Cell-cycle phase per cluster/timepoint | Seurat `CellCycleScoring`; plot via `plot_cluster_summary.py` | [orig] |
| S8 A | **Cluster-17 de-differentiation program scores** (heatmap: EMT peak; CSC/neural-crest/MSC elevated but peak in neighbors; cl17-vs-rest MWU/BH stars) *(redesign of orig Supp Fig 3)* | `scripts/python/plot_cl17_stemness_programs.py` (+ `data/cl17_stemness_signatures.tsv`) | [HAVE] |
| S8 B | **Cluster-17 curated marker dot plot** (de-duplicated, separate TF block; ★ = BH-significant positive DEG, read enrichment from ★ not color) *(redesign)* | `scripts/python/plot_marker_dotplot.py --markers data/cl17_stemness_markers.tsv --deg-table data/supp_table3_findallmarkers.tsv` | [HAVE] |
| S9 | **Cluster-17 TP63 + basal signature** | `plot_cl17_basal_vs_emt.py` (+ `data/cl17_basal_emt_signatures.tsv`) | [HAVE] |
| S10 | GSEA clusters 4/6 vs cl0 and vs cl17 | `per_cluster_gsea.py` (+ `cluster4_emt_gradient.py`) | ~ |
| S11 | Smoothed single-cell expression heatmap of 50 lineage-associated TFs along velocity-pseudotime (`data/s11_panelA_tfs.tsv`) | `plot_pseudotime_gene_heatmap.py` | [HAVE] |
| S12 | **Trajectory tracing scores + temporal stats** | `jonckheere_trend.py` + `plot_pseudotime_qc.py` + `extract_trajectory_summary.py` | [HAVE] |
| S13 | RSS regulon map **+ canonical NEPC TFs** | `plot_scenic_regulon_map.py` + `scenic_top_regulons.py` | [HAVE] |
| S14 | **SCENIC↔decoupleR consensus** | `scenic_decoupler_consensus.py` | [HAVE] |
| S15 | **Chromatin-surrogate evidence** — (A) full 17-TF Shrestha ATAC footprint atlas (companion to curated 5D), (B) Baca HOMER motif enrichment (NEPC- vs PRAD-specific regions; NE + neural-crest tier), (C) ChEA3 integrated TF–target rank map. decoupleR dropped (=5B/S14) | `plot_shrestha_chromatin.py` (full, no `--curated`) + `chromatin_surrogate/` HOMER pipeline (`data/chromatin_surrogate/public/figures/fig5C_chromatin_motif` → **S15B**) + ChEA3 `plot_chea3_map.py` (`fig5_S15C_chea3_map` → **S15C**) | [HAVE] |
| S16 | **GEP correspondence** across 18 clusters + 8 timepoints (anchors GEP bulk-validation S17/S18) | `plot_gep_correspondence.py` | [HAVE] |
| S17 | GEP7 in dormant PDX + CRPC/NEPC | cNMF (orig) | [orig] |
| S18 | GEP1 (ASCL1−) / GEP6 (ASCL1+) + AR markers | cNMF (orig) | [orig] |
| S19 | Gao/Dong dataset processing | dataset processing (orig); markers via `plot_marker_dotplot.py` | [orig]/~ |
| S20 | Li/Wang dataset processing | dataset processing (orig) | [orig] |
| S21 | **Integration QC** (iLISI/cLISI/PERMANOVA) | `integration_lisi_kbet.py` + `integration_permanova_knn.py` | [HAVE] |
| S22 | **Clustering-resolution justification** (0.55) **+ cl17 directed retrieval** | `harmony_cluster_agreement.py` (+ `compare_resolution_trajectory.py`) | [HAVE] |
| S23 | **Patient-level signature transfer** | `make_signatures.py` → `clinical_signature_validation.py` | [HAVE] |
| S24 | **9- vs 10-patient (Drop-seq P3) integration sensitivity** | `integration_lisi_kbet.py` (9-pt) + `cluster_origin_correspondence.py` | [HAVE] |
| S25 | **Cluster-4 signature vs progression-free interval, stratified by ERG fusion** (TCGA-PRAD) | `tcga_cl4_erg_bcr.py` | [HAVE] |
| S26 | **Truncal deletion-type TMPRSS2:ERG fusion (WES): present pre-castration, retained through NEPC** | `plot_tmprss2_erg_deletion.py` | [HAVE] |

---

## Supplementary tables (canonical T1–T9)

| Table | Content | Source / Script | Status |
|-------|---------|-----------------|--------|
| T1 | Sample / timepoint metadata | `annotate_timepoint.py` / metadata | [orig] |
| T2 | Cell-type & lineage marker gene lists | `data/Supplementary_Table2.txt`; `data/fig2c_*.tsv` | [HAVE] |
| T3 | Per-cluster DEGs (FindAllMarkers) | Seurat | [orig] |
| T4 | SCENIC: 200 regulons + RSS | `scenic_top_regulons.py` (pySCENIC) | [HAVE] |
| T5 | Doublet summary **+ per-cluster rate + Fisher** | `doublet_robustness.py` | [HAVE] |
| T6 | **KRAS_UP / KRAS_DN leading-edge genes** | `geneset_leading_edge.py` (+ `kras_lineage_overlap.py`) | [HAVE] |
| T7 | **NED→NE enrichment per cohort** | `annotate_validation_phenotypes.py --report` / `plot_ned_enrichment.py` | [HAVE] |
| T8 | **Integration QC metrics** | `integration_lisi_kbet.py` + `integration_permanova_knn.py` + `harmony_cluster_agreement.py` | [HAVE] |
| T9 | **Methods parameter table** | `supplementary/SUPPLEMENTARY_TABLE_T9_parameters.md` (+ `config/params.yaml`) | [HAVE] |

---
