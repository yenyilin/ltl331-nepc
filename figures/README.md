# Figures

Rendered figure panels and their source data. Everything here is **generated** —
regenerate from `../scripts/` (Python) and `../R/` following `../RUN.md`; the
script → panel → path map is in `../FIGURES.md`.

## Layout

```
figures/
├── main/          # main-text figure panels + assembled composites
├── supp/          # supplementary figures (S1–S25)
└── source_data/   # the numeric data behind each panel
```

## Naming convention

Name every output by its **figure number** so a panel is findable at a glance.

| kind | pattern | example |
|---|---|---|
| main panel | `fig<N><letter>_<slug>.{pdf,png}` | `fig6B_metaneighbor.pdf` |
| main composite | `Fig<N>.pdf` | `Fig6.pdf` |
| supplementary | `S<N>_<slug>.{pdf,png}` | `S6_wes_infercnv_concordance.pdf` |
| source data | `S<N>_<slug>.tsv` / `fig<N>_<slug>.tsv` | `source_data/S6_per_chrom_correlations.tsv` |

- Keep **both `pdf` and `png`**: PDF is the vector version for the manuscript;
  PNG renders inline on GitHub.
- `supp/` is flat (supplementary figures are mostly single-panel). For a
  many-panel *main* figure, an optional `main/fig<N>/` subfolder may hold the
  individual panels, with the composite as `main/Fig<N>.pdf`.
- `source_data/` mirrors the figure numbers — one table per panel that has
  underlying numbers (enrichment stats, correlations, scores, …).

## Rule

Do not hand-edit files here. If a panel needs changing, change the script and
re-run — the figure is an artifact, the script is the source of truth.
