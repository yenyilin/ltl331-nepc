# Data Availability

---

## Processed objects (GEO)

| Object | Format | DOI / link |
|--------|--------|-----------|
| LTL331/R Seurat object | GSE297328_LTL331.RDS |[GSE297328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297328) |
| LTL331/R h5ad object | GSE297328_ltl331_annotated.h5ad |[GSE297328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297328) |
| LTL331/R velocity h5ad object | GSE297328_ltl331_velocity.h5ad |[GSE297328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297328) |
| LTL331/R and patient integrated dataset | GSE297328_ltl331_harmony_integrated.h5ad | [GSE297328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297328) |

## Raw sequencing (GEO / NCBI)

| Data | Accession | Access |
|------|-----------|--------|
| LTL331/R scRNA-seq | [GSE297328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297328) | publicly available |

## Public datasets reused

| Dataset | Accession |
|---------|-----------|
| Gao et al. (clinical scRNA-seq) |[GSE137829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137829) |
| Li et al. (clinical scRNA-seq) | GSA-Human HRA002145 |
| Bulk RNA-seq compendium (316 CRPC and 19 NEPC samples used) | Processed, batch-corrected compendium of [Bolis et al. 2021](https://doi.org/10.1038/s41467-021-26840-5). Source accessions of the samples used: dbGaP [phs000915](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000915), [phs000673](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000673), [phs000909](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000909) (controlled access); GEO [GSE126078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126078), [GSE118435](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118435); ENA [PRJEB21092](https://www.ebi.ac.uk/ena/browser/view/PRJEB21092); TCGA-PRAD |
| TCGA-PRAD (PanCancer Atlas; cluster-4 × ERG recurrence analysis) | cBioPortal study [prad_tcga_pan_can_atlas_2018](https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018) |
| NEPC/PRAD LuCaP PDX ChIP-seq (HOMER motif analysis) | GEO [GSE161948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161948) ([Baca et al. 2021](https://doi.org/10.1038/s41467-021-22139-7)) |
| Metastatic CRPC ATAC-seq atlas (TF footprint occupancy) | Processed supplementary tables of [Shrestha et al. 2024](https://doi.org/10.1158/0008-5472.CAN-24-0890) (West Coast Dream Team cohort, n = 70) |
| Normal prostate reference (inferCNV) | d17 normal adult prostate 3 ([Henry et al. 2018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120716)) |

## Code

- This repository: `https://github.com/yenyilin/ltl331-nepc` — see `FIGURES.md`
  for the script behind each panel and `docs/methods_parameter_table.md` for
  all parameters.
- Archived on Zenodo: https://doi.org/10.5281/zenodo.23159059 (concept DOI,
  always resolves to the latest version; v1.0.0 = 10.5281/zenodo.23159060).
---
