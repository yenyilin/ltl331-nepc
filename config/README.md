# Config

`params.yaml` is the **single source of truth** for every threshold, seed,
and pinned software version used in the analysis. No script should hardcode a
number that belongs here.

- **Human-readable mirror:** `../docs/methods_parameter_table.md` is derived
  from this file. Update that table accordingly since they're not auto-synced.
- **Usage:** scripts take the equivalent values as **CLI flags** with the
  `params.yaml` value as the documented default (check each script's `--help`)
  rather than reading the YAML directly.

## Sections

| Key | Covers |
|-----|--------|
| `versions` | pinned tool versions used across the study (Cell Ranger, Seurat, UCell, inferCNV, fgsea, scVelo, CellRank, pySCENIC, cNMF, Harmony, sva; R + Python) |
| `qc` | PDX species filter, mito %, MAD/UMI/gene-count cutoffs |
| `clustering` | PCs, cluster count, seed |
| `infercnv` | reference, subsampling, cutoff, HMM settings |
| `trajectory` | velocity-graph neighbors/PCs, PAGA threshold, scVelo mode |
| `cellrank` | GPCCA estimator, kernel weights, terminal-state definitions (original vs revision) |
| `scenic` | subsampling, iterations, regulon-count thresholds |
| `cnmf` | GEP factorization range/iterations |
| `public_integration` / `joint_integration_revision` | cross-cohort Harmony settings (patient-level vs revision joint integration) |
| `bulk_rnaseq` | CRPC/NEPC cohort size, batch correction method |
| `seeds` | global + Leiden seeds |
| `data_objects` | logical paths/keys for the three canonical AnnData objects — full contract in `../data/objects/MANIFEST.md` |

## Before submission

The whole point of this file is that it's traceable to what was really run. 
Then re-derive `../docs/methods_parameter_table.md` to match.
