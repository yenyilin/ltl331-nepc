# Data

No data files are stored in this repository (single-cell objects are GBs; GitHub
caps at 100 MB/file). All data is hosted externally — see `../DATA_AVAILABILITY.md`
for accessions and links.

## Expected local layout after fetching (gitignored)

```
data/
├── objects/          # the 3 canonical AnnData inputs — see objects/MANIFEST.md
│   ├── ltl331_annotated.h5ad           # (1) 51,726 cells; .X scaled, lognorm in .raw
│   ├── ltl331_velocity.h5ad            # (2) velocity layers + CellRank2 fate obs
│   ├── ltl331_harmony_integrated.h5ad  # (3) LTL331+Gao+Li joint Harmony
│   ├── MANIFEST.md                # per-object schema contract (committed)
│   └── CHECKSUMS.txt              # sha256 of each object (generate from the GEO objects, then commit)
├── tables/           # small committed inputs (regenerate or version-control)
│   ├── genesets/     # e.g. hallmark_v7.5.1.tsv (from MSigDB GMT; see config)
│   ├── cnv/          # nexus_*, infercnv_*, WES CNV calls (TSV) — see docs/nexus_cnv_figure.md
│   └── scenic/       # regulon / RSS tables (TSV)
├── doubletfinder/    # Fig S1 before/after count matrices — see doubletfinder/README.md
├── raw/              # GEO GSE297328 (CellRanger outputs) if re-running from scratch
└── external/         # Gao (GSE137829), Li (HRA002145), bulk (Bolis et al. 2021 compendium; see DATA_AVAILABILITY.md),
                      # inferCNV normal reference, SCENIC motif databases
```

`objects/` holds large gitignored `.h5ad`; only its `MANIFEST.md` is committed
(add `CHECKSUMS.txt`, generated from the deposited GEO objects, to let a fetch be
verified). `tables/` holds small derived inputs that *can* live in the repo (so
a figure can be remade without the GB-scale objects) — commit them or document how
to regenerate.

## The three objects (why three)

Each analysis family binds to one object; they have **different internal schemas**
(matrix normalization, layers, obs keys). The mapping and the load rules are in
`objects/MANIFEST.md`; logical paths/keys are mirrored in
`../config/params.yaml › data_objects`.

| object | drives | load note |
|--------|--------|-----------|
| `ltl331_annotated.h5ad` | Fig 2/3, markers, GSEA, SCENIC, decoupleR | score from `.raw` (`.X` is scaled) |
| `ltl331_velocity.h5ad` | Fig 4 trajectory / fate | velocity layers + CellRank obs |
| `ltl331_harmony_integrated.h5ad` | Fig 6 integration / clinical | check `.X` min before pseudobulk |

## Fetch (examples)

```
# processed objects from GEO GSE297328 (public):
# wget -O data/objects/ltl331_annotated.h5ad \
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297328/suppl/GSE297328_ltl331_annotated.h5ad
# wget -O data/objects/ltl331_velocity.h5ad \
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297328/suppl/GSE297328_ltl331_velocity.h5ad
# wget -O data/objects/ltl331_harmony_integrated.h5ad \
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297328/suppl/GSE297328_ltl331_harmony_integrated.h5ad
# sha256sum -c data/objects/CHECKSUMS.txt        # verify integrity (once CHECKSUMS.txt is generated — see below)

# scripts read paths via --h5ad / CLI args or config; do not hardcode.
```

Scripts take data paths as arguments (e.g. `--h5ad`), so point them at wherever you
place these files.
