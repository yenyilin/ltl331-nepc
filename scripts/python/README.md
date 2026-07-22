# Python scripts

The revision analysis code, everything that produces a main-text or
supplementary panel, plus the compute steps those panels depend on. Flat
directory (no subpackages) except `chromatin_surrogate/`, which is a
standalone sub-pipeline with its own environment. This README groups the 
same scripts by role so you can find the one you want.

**Don't move or rename files here.** Every script is launched as
`python scripts/python/foo.py` and imports siblings (`from _utils import ...`)
via `sys.path[0]` — there's no package/`__init__.py`. Splitting this into
subfolders would break those imports.

For "which script makes panel X" see `../../FIGURES.md`. For "in what order
do I run these" see `../../RUN.md`.

## Shared helper

- `_utils.py` — plot style (`setup_style()`: embedded TrueType fonts so PDF
  text stays editable; `image.interpolation='nearest'` so thin color-strip
  panels don't render blank in the PDF backend), `savefig()`, and the sort
  keys (`numkey`, `timepoint_key`) used across most plotters. Import this
  instead of re-implementing any of the above.
- `make_cluster_colors.py` — derives the 18-shade phenotype-banded cluster
  palette every figure uses.
- `assemble_figure.py` — composes label-free single-panel PDFs into the final
  multi-panel figure mosaics.

## By role

**Doublet robustness** (Supp S1, R1#1)
`doublet_unify.py` → `doublet_robustness.py`

**CNV: WES / inferCNV** (Supp S2/S4/S5/S6, R1#1/R2#2)
`nexus_cnv_loader.py`, `infercnv_loader.py`, `plot_nexus_cnv_curves.py`,
`plot_nexus_cnv_heatmap.py`, `compare_infercnv_vs_wes.py`,
`infercnv_clonal_check.py`, `cnv_driver_locus_tracks.py`
— see `../../docs/nexus_cnv_figure.md` for the Nexus curve figure specifically.

**Integration + clinical validation** (Fig 6, Supp S21–S25, R1#3/R2#7)
`harmony_cluster_agreement.py`, `compare_resolution_trajectory.py`,
`plot_origin_vs_harmony.py`, `cluster_origin_correspondence.py`,
`patient_level_stats.py`, `integration_lisi_kbet.py`,
`integration_permanova_knn.py`, `annotate_validation_phenotypes.py`,
`ned_ne_enrichment.py`, `plot_ned_enrichment.py`, `terminal_metaneighbor.py`,
`plot_terminal_metaneighbor.py`, `terminal_reference_project.py`,
`plot_terminal_projection.py`, `plot_fig6a_label_transfer.py`,
`make_signatures.py`, `clinical_signature_validation.py`

**Trajectory: CellRank / velocity / PAGA** (Fig 4, R1#5/R2#6)
`cellrank_bifurcation.py`, `paga_trajectory.py`, `velocity_dynamical.py`,
`plot_terminal_states.py`, `plot_fate_umap.py`, `cluster17_absorption.py`,
`cluster17_transition_flux.py`, `jonckheere_trend.py`,
`verify_paga_obligate.py`, `plot_pseudotime_qc.py`

**Regulon / GEP / cell-state**
`decoupler_tf_activity.py`, `plot_gep_correspondence.py`,
`plot_phase_composition.py`, `plot_marker_dotplot.py`

**Fig 5: SCENIC / decoupleR / regulon maps** (R1#4, R2#1)
`plot_scenic_regulon_map.py`, `plot_regulon_driver_map.py`,
`plot_tf_activity_heatmap.py`, `plot_branch_tf_tracks.py`,
`plot_shrestha_chromatin.py`, `plot_onecut2_expression.py`,
`plot_chea3_map.py` — plus the standalone `chromatin_surrogate/` sub-pipeline
(own README/env; produces the chromatin anchor these scripts render).

**Fig 2: UMAP / markers / timepoints**
`plot_umap_panels.py`, `plot_enrichment_violin.py`,
`plot_timepoint_panel.py`, `plot_cl17_basal_vs_emt.py`

**Fig 3: proportions / DEG / GSEA / pathways** (R2#5)
`plot_cluster_proportions.py`, `plot_fig3b_timepoint_ridge.py`,
`plot_fig3c_deg_heatmap.py`, `fig3d_nes_heatmap.py`, `per_cluster_gsea.py`,
`plot_fig3e_cl17_pathways.py`, `plot_fig3e_cl4_12_dynamics.py`,
`plot_fig3f_proliferation_violin.py`, `pathway_revision.py`,
`kras_lineage_overlap.py`, `make_genesets.sh` (builds the pinned MSigDB
v7.5.1 gene-set TSVs these consume — the one shell script in this directory)

**Cluster-17 identity / stemness / novelty** (Supp S8, tuft-negativity panel)
`plot_cl17_stemness_programs.py`, `plot_pou2f3_tuft_panel.py`

## What's deliberately not here

Notably absent on purpose: the superseded Fig 6A pair (Harmony co-embedding design, 
replaced by the embedding-free panels above) and any script still under active debugging.
