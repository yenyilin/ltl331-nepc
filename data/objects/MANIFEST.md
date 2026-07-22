# Canonical data objects — schema manifest

The analysis runs off **three** AnnData objects. Each has a *different internal
schema*; scripts must read the right matrix and keys per object. This manifest is
the contract. (Objects are gitignored and hosted on figshare — see
`../../DATA_AVAILABILITY.md`. Logical paths and key names are mirrored in
`config/params.yaml › data_objects`.)

> **Load rule of thumb:** never assume `.X` is lognorm. Check `adata.X.min()`
> (negative ⇒ scaled) and prefer `.raw` / a named lognorm layer for any gene-set
> scoring or mean-expression. Object (1) is confirmed scaled-`.X`.

---

## (1) `ltl331_annotated.h5ad` — base expression object

51,726 cells. Drives clustering, composition, markers/DEG, GSEA/pathway scoring,
proliferation, SCENIC input, decoupleR.

| slot | content | use |
|------|---------|-----|
| `.X` | **SCALED / z-scored**, HVG subset (min ≈ −6.93, max clipped 10.0 = Seurat `ScaleData`) | **do NOT score on this** |
| `.raw` | **lognorm, full 26,175 genes** | score from here (`--use-raw` / `use_raw=True`) |
| `.layers` | *(empty)* | — |
| `.obs` | `seurat_clusters` (int 0–17 = manuscript IDs), `week_num` (0/2/4/8/12/16/20/22), `timepoint` (`preCx`,`wk2`…`wk22R` — note `R` suffix), `is_regrowth`, `phase_castration`, `orig.ident`, `pretty`, `cellid`, `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `*_UCell` signature scores | cluster key = `seurat_clusters`; timepoint key = `week_num` |

**Consumers:** `pathway_revision.py` (`--use-raw`), `tier0_verification.py`,
`decoupler_tf_activity.py`, `cluster4_emt_gradient.py`, `cluster6_split.py`,
`plot_fig3c_deg_heatmap.py`, `plot_fig3f_proliferation_violin.py`,
`plot_cluster_summary.py`, `extract_per_timepoint_profile.py`, `make_signatures.py`.

---

## (2) `ltl331_velocity.h5ad` — velocity + CellRank2 object

Same cells, carrying RNA-velocity layers and CellRank2 fate results. Drives PAGA,
velocity, pseudotime, fate/absorption (Fig 4), trajectory schematic, pt-QC.

| slot | content | use |
|------|---------|-----|
| `.layers` | `spliced`, `unspliced`, `velocity`, `Ms`, `Mu` | velocity kernel inputs |
| `.obs` | `seurat_clusters`, `week_num`, `velocity_pseudotime`, `latent_time`, CellRank macrostate + absorption-probability columns `<confirm exact names>` (e.g. `P_ASCL1pos`,`P_ASCL1neg`,`macrostates`) | cluster key = `seurat_clusters`; pt = `velocity_pseudotime` |
| `.obsm` | `X_umap`, velocity embedding `<confirm>` | plotting |
| `.uns` | PAGA, velocity graph, CellRank estimator `<confirm>` | trajectory |
| `.X` / `.raw` | `<confirm — likely shares base lognorm in .raw>` | — |

**Consumers:** `cellrank_bifurcation.py`, `extract_trajectory_summary.py`,
`plot_pseudotime_qc.py`, `plot_trajectory_schematic.py`,
`plot_cluster_temporal_dynamics.py`, `compare_resolution_trajectory.py`.

> `<confirm>` items: fill exact obs/obsm/uns key names when the object is deposited
> (or inspect with `adata.obs.columns`, `adata.layers.keys()`, `adata.obsm.keys()`).

---

## (3) `harmony.annotated.h5ad` — joint Harmony integration

LTL331 + Gao + Li, jointly embedded (Leiden 0.55). Drives integration QC,
cross-cohort correspondence, clinical co-mapping (Fig 6), patient-level stats.

| slot | content | use |
|------|---------|-----|
| `.obsm` | `X_pca_harmony` (integrated embedding), `X_umap` | LISI/kBET, co-mapping |
| `.obs` | cohort/batch key `<confirm: cohort/dataset>`, patient key `<confirm: patient/sample>`, joint cluster `leiden_0.55`, original LTL331 label `seurat_clusters` carried through (NA for clinical cells), transferred-label column `<confirm>` | batch = cohort; patient = `<confirm>` |
| `.X` / `.raw` | **CHECK normalization** — `patient_level_stats.py` reads `.X` unless `--layer` given and has **no `.raw` fallback**; if `.X` here is scaled, pass `--layer <lognorm>` | pseudobulk signatures |

**Consumers:** `harmony_cluster_agreement.py`, `cluster_origin_correspondence.py`,
`plot_origin_vs_harmony.py`, `patient_level_stats.py` (**verify matrix**),
`transfer_dynamics_to_clinical.py`, `clinical_bifurcation_balance.py`.

> Before running `patient_level_stats.py`: `print(adata.X.min())` on this object —
> if negative, supply `--layer <lognorm-layer>`.

---

## Integrity

Record checksums + provenance when depositing (so a fetched object is verifiably
the one the figures were made from):

| object | figshare DOI | sha256 | built by | date |
|--------|--------------|--------|----------|------|
| ltl331_base.h5ad | `<TBD>` | `<TBD>` | `<pipeline/commit>` | `<TBD>` |
| ltl331_velocity_cr2.h5ad | `<TBD>` | `<TBD>` | `<commit>` | `<TBD>` |
| ltl331_harmony.h5ad | `<TBD>` | `<TBD>` | `<commit>` | `<TBD>` |

`sha256sum data/objects/*.h5ad > data/objects/CHECKSUMS.txt`
