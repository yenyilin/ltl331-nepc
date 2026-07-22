# Supplementary Tables — data files

Machine-readable supplementary tables for the LTL331/R NEtD revision. Each `.tsv`
carries a `#`-commented provenance header. Full legends are in the manuscript
supplement (`response/SUPPLEMENTARY_clean.md`, "Supplementary Tables" section).

| File | Supplementary Table | Reviewer | Source data | Assembly script |
|------|---------------------|----------|-------------|-----------------|
| `T1.tsv` | **T1.** Sample / timepoint metadata | — | `data/doubletfinder/unified/doublet_by_timepoint.tsv` | — |
| `T2.tsv` | **T2.** Cell-type / lineage marker gene lists (canonical + Fig 2D panel + Fig 3 cluster markers) | R2#4 | `data/Supplementary_Table2.txt`, `data/fig2c_main.tsv`, `data/fig3_datadriven.tsv` | — |
| `T3.tsv` | **T3.** Per-cluster DEGs (FindAllMarkers) | — | `data/supp_table3_findallmarkers.tsv` | Seurat `FindAllMarkers` |
| `T4.tsv` | **T4.** SCENIC 200 regulons + RSS | R1#4 | `data/scenic/8samples_exp_15_all_all_cluster_RSS.tsv` | pySCENIC (20-run consensus) |
| `T5.tsv` | **T5.** Per-cluster doublet-removal summary + Fisher enrichment | R1#1 | `data/doubletfinder/unified/doublet_by_cluster.tsv` + `cluster_enrichment.tsv` | `doublet_unify.py` → `doublet_robustness.py` |
| `T6.tsv` | **T6.** HALLMARK_KRAS_SIGNALING_UP/_DN leading-edge genes per NE cluster (+ `NE_lineage_marker`) | R2#5 | `data/kras_lineage_overlap/kras_lineage_overlap_SuppTableX.tsv` | `scripts/kras_lineage_overlap.py` |
| `T7.tsv` | **T7.** NED→NE epithelium enrichment per cohort (2×2 counts, Fisher OR, p) | R1#3 | `data/integration_qc/validation_ne_enrichment.tsv` | `scripts/integration_permanova_knn.py` |
| `T8.tsv` | **T8.** Integration QC metrics (iLISI/cLISI, PERMANOVA, dispersion, cl17 directed-kNN, resolution sweep) | R1#3 | `data/integration_qc/*` + `data/agreement/metrics_summary.tsv` | `integration_permanova_knn.py` + `harmony_cluster_agreement.py` |

**T9** (methods parameter table) is auto-generated from `config/params.yaml` by the
assembler below (78 params: versions, QC, clustering, inferCNV, trajectory, CellRank).

## Assembling the submission workbook

`python3 scripts/assemble_supp_tables.py` combines T1–T8 (+ generated T9) into
**`data/supp_tables/supplementary_tables.xlsx`**, a Contents sheet plus one sheet per
table (title + legend + provenance rows, then the data).  
## Notes

- **T1/T5** counts are barcode-harmonized across the pre/post-doublet objects
  (`data/doubletfinder/unified/`); totals 51,726 → 49,568 analyzed cells (rate 4.17%).
  The older per-cluster/-timepoint rate files are deprecated (they produced
  impossible negative doublet counts). Cluster 17 is depleted (OR 0.36).
- **T2** is long format (`source, group, gene`); `source` distinguishes the canonical
  annotation panels from the Fig 2D program panel and the Fig 3 per-cluster markers.
- **T3** is one row per (cluster, gene); Bonferroni-adjusted `p_val_adj`.
- **T4** is long format (`cluster, regulon, RSS, rss_rank_in_cluster`); regulon (+)/(−)
  = activating/repressing; rank 1 = most cluster-specific.
- **T7** 95% CIs for the odds ratios are in main Fig. 6D (Gao 60–300; Li 24–41). The
  Li Fisher p underflows to 0 (< 1×10⁻³⁰⁰).
- **T8** is tidy long (`section, metric, condition, value`); `naive_PCA` = before Harmony,
  `Harmony` = after. Selected clustering resolution = 0.55.

## Related repo data files (former original Supp Tables 6/7/8)

- `data/cellrank2_lineage_drivers.tsv` — CellRank2 lineage-driver correlations (backs Supp Fig S11)
- `data/gep_top100.tsv` — top-100 genes per gene-expression program (pairs with main Table 1)
- `data/cluster_timepoint_cell_counts.tsv` — analyzed cell counts per cluster × timepoint (backs Fig 3A)
