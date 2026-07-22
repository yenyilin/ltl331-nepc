# Nexus CNV figure (`plot_nexus_cnv_curves.py`)

Curve-style copy-number visualisation of the LTL331/R Nexus WGS cohort — the
cohort gain/loss frequency plot plus per-sample CNV profiles, with driver-gene
annotations and an LOH / allelic-imbalance track. Build **GRCh38 / hg38**.

## Input data

All coordinates are **hg38**. Sample names are the public (sanitized) IDs — the
normalization-reference suffix (`vs …`) is not stored.

| file | contents |
|---|---|
| `data/nexus_calls.tsv` | CN segments — `sample`, `Chromosome`, `Start`, `End`, `Value` (−2…+2), `value_label` |
| `data/nexus_zygosity.tsv` | LOH / allelic-imbalance segments (same schema; `value_label` = `imbalance` / `hom_high` / `hom_low`) |
| `data/nexus_descriptors.tsv` | per-sample metadata — `sample`, `timepoint`, `series` |
| `data/nexus_sample_order.tsv` | row order + display labels — `sample`, `display_name`, `order` |
| `data/driver_genes_hg38.tsv` | loci to annotate — `gene`, `chrom`, `start`, `end` |

These TSVs are regenerated from the raw Nexus export by
`scripts/nexus_cnv_loader.py`, which writes public-safe sample names by default.

## Usage

```bash
python3 scripts/plot_nexus_cnv_curves.py \
  --calls-tsv       data/nexus_calls.tsv \
  --descriptors-tsv data/nexus_descriptors.tsv \
  --bin-size        1000000 \
  --mode            combined \
  --build           hg38 \
  --genes           data/driver_genes_hg38.tsv \
  --zygosity        data/nexus_zygosity.tsv \
  --sample-order    data/nexus_sample_order.tsv \
  --outdir          nexus_plots \
  --plot-format     pdf png
```

→ `nexus_plots/cnv_curves_combined.{pdf,png}`

## Options

| flag | default | description |
|---|---|---|
| `--calls-tsv` | *(required)* | CNV segment calls |
| `--descriptors-tsv` | – | sample metadata (timepoint-aware sorting) |
| `--bin-size` | `1000000` | genome bin width in bp (keep 1 Mbp to match the concordance figure) |
| `--mode` | `per_sample` | `per_sample` · `frequency` · `combined` |
| `--build` | `hg38` | chromosome sizes + x-axis label (`hg38` or `hg19`) |
| `--genes` | – | driver-gene TSV → dashed markers + gene labels down every row |
| `--zygosity` | – | LOH/imbalance TSV → coloured band beneath each per-sample row (purple = allelic imbalance, orange-brown = LOH) |
| `--sample-order` | – | row order + biological display labels |
| `--filter-samples` | – | restrict to a subset of samples |
| `--title` | – | figure suptitle |
| `--save-matrix` | off | also write the painted sample × bin integer matrix |
| `--plot-format` | `pdf` | `pdf` / `png` / `svg` (space-separated) |

## Modes

- **`combined`** — cohort gain/loss frequency mirror plot on top, per-sample CNV
  profiles below, shared genome x-axis. Recommended for the supplement.
- **`per_sample`** — per-sample small multiples only.
- **`frequency`** — cohort recurrence mirror plot only.

Driver-gene markers and the LOH/imbalance band render in `combined` and
`per_sample`; `frequency` shows the gene markers only.

## Outputs

| mode | file |
|---|---|
| `combined` | `nexus_plots/cnv_curves_combined.{pdf,png}` |
| `per_sample` | `nexus_plots/cnv_curves_per_sample.{pdf,png}` |
| `frequency` | `nexus_plots/cnv_curves_frequency.{pdf,png}` |
| `--save-matrix` | `nexus_plots/cnv_curves_matrix.tsv` |

## Example variants

```bash
# cohort frequency mirror plot with driver genes
… --mode frequency --genes data/driver_genes_hg38.tsv

# per-sample profiles with LOH band + biological row order
… --mode per_sample --zygosity data/nexus_zygosity.tsv --sample-order data/nexus_sample_order.tsv

# the four single-cell–profiled timepoints only
… --filter-samples LTL331_preCX1 LTL331_16wk1 LTL331_20wk LTL331_22wk2
```
