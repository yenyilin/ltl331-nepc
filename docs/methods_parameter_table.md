# Methods Parameter Table

Human-readable transparency table. Values pre-filled from the manuscript Methods;
Machine-readable source: `config/params.yaml`.

> Capture versions/seeds from the machine that produced the results
> (`uv pip freeze`, scanpy `print_versions()`) rather than transcribing —
> see `env/`.

---

## Software versions

| Tool | Version | Step |
|------|---------|------|
| Cell Ranger | 6.0.0 | alignment, aggr |
| Seurat | 4.1.1 (Gao/Li: 4.3.0) | QC, clustering, DEG |
| UCell | 2.6.2 | signature scoring |
| inferCNV | 1.15.0 | CNV |
| msigdbr / fgsea | 7.5.1 / 1.30.0 | GSEA |
| scVelo | 0.3.4 | RNA velocity |
| CellRank | 2.0.7 | fate probabilities |
| pySCENIC | 0.12.0 | regulons |
| cNMF | 1.4.1 | GEPs |
| Harmony | 1.1.0 | integration |
| sva (ComBat) | 3.50.0 | bulk batch correction |
| R / Python | 4.2.1 / 3.11.9 | all |

## QC and filtering

| Parameter | Value |
|-----------|-------|
| Human read fraction (PDX species filter) | > 0.80 |
| Mitochondrial transcript % | < 25% |
| UMI / gene counts | within 3 MAD of median |
| Minimum UMIs per cell | ≥ 1,000 |
| Min cells per gene | > 10 |
| Doublet removal retention | > 95% (Supp Table 9) |

## Clustering / dimensionality reduction

| Parameter | Value |
|-----------|-------|
| PCs retained | 50 (knee point) |
| Clusters (LTL331) | 18 |
| Normalization | log-normalization (Seurat defaults) |
| HVGs | Seurat defaults / `<TBD>` |
| Random seed | `<TBD — record actual>` |

## inferCNV

| Parameter | Value |
|-----------|-------|
| Reference | luminal epithelial, d17 normal adult prostate 3 |
| Subsample per group | 2,000 cells |
| cutoff | 0.1 (10x recommended) |
| hclust_method | ward.D2 |
| noise_filter | 0.01 |
| leiden_resolution | 0.001 |
| HMM | TRUE, type i6 |
| Groupings | by sample and by cluster; cluster_by_groups TRUE/FALSE |

## Trajectory / velocity / fate

| Parameter | Value |
|-----------|-------|
| PAGA edge weight threshold | 0.09 (Fig 4A) |
| scVelo genes | top 2,000; min 20 counts spliced & unspliced |
| scVelo mode | dynamical |
| velocity pseudotime | `tl.velocity_pseudotime` |
| CellRank estimator | GPCCA |
| CellRank kernel (revision) | 0.8·VelocityKernel + 0.2·ConnectivityKernel |
| CellRank macrostates (revision) | 8 terminal states recovered `{0,1,4_1,4_2,6,7,9,10}`; 4 NE fates (clusters 1, 7, 9, 10) |
| Terminal states (original) | cluster 7 (sole NE terminal) |
| Terminal states (R1#5 revision) | ASCL1+ (cluster 10) & ASCL1− (cluster 7) |
| TF annotation source | cisTarget TF list (fetched Aug 2022) |

## SCENIC

| Parameter | Value |
|-----------|-------|
| Subsampled cells | 7,878 (≤500 per cluster) |
| Iterations | 20 |
| TF kept if reported | ≥ 16 / 20 runs |
| Target kept if in | ≥ 80% of runs identifying that TF |
| Active regulons (total) | 200 |
| Motif database | Motif2TF V9 |

## cNMF

| Parameter | Value |
|-----------|-------|
| k range scanned | 5–20 |
| iterations | 20 |
| HVGs | 2,000 |
| density threshold | 0.03 |
| GEPs selected | 7 |
| Top genes per GEP | 100 |

## Public dataset integration (Gao / Li)

| Parameter | Value |
|-----------|-------|
| Within-cohort integration | Harmony 1.1.0 |
| Epithelial gating (Gao) | author annotation |
| Epithelial gating (Li) | marker genes (Supp Table 2) |
| Signature scoring | UCell, top-20 DEG (log2FC>0.5, FDR<0.05) |
| Significance | Wilcoxon rank-sum (cluster vs rest) |

## Joint integration 

| Parameter | Value |
|-----------|-------|
| Method | Harmony (primary) |
| Clustering resolution | Leiden 0.55 → 27 clusters |
| Harmony theta / max_iter | `<TBD — record>` |
| Random seed (leiden) | 0 |

## Bulk RNA-seq

| Parameter | Value |
|-----------|-------|
| Cohort | 316 CRPC + 19 NEPC (E-MTAB-9930) |
| Batch correction | ComBat (sva 3.50.0), PolyA+ reference batch |
| GEP scoring | top-30 genes per GEP, hierarchical clustering |

---
