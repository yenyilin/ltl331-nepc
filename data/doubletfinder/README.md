# Doublet sensitivity inputs (Supplementary Fig. S1, R1 #3)

Two harmonized per-(timepoint × cluster) cell-count matrices, **identical layout**
(`timepoint` rows × `cluster_0 … cluster_17` columns):

- `counts_original.tsv` — cells **before** DoubletFinder
- `counts_valid.tsv`    — the same grid **after** doublet removal

## What the step shows
DoubletFinder is run per timepoint on the full object; these two matrices are the
per-(timepoint, cluster) cell counts before and after removal. Differencing them
(`doublets = original − valid`) gives a per-cluster / per-timepoint doublet rate,
which is what Supplementary Fig. S1 uses to show the load-bearing states, in
particular the EMT/bridge cluster 17, are **not** doublet artifacts.

## Flow
```
counts_original.tsv ┐
                    ├─▶ doublet_unify.py ──▶ unified/ (matrices, long table,
counts_valid.tsv    ┘        │                per-cluster, per-timepoint, flags)
                             ▼
                    doublet_robustness.py ──▶ doublet_robustness.{pdf,png}  (= Fig. S1)
                                              + doublet_robustness_stats.tsv
                                              + cluster_enrichment.tsv (Fisher OR + BH-q)
```

## Regenerate
```
python scripts/python/doublet_unify.py --dir data/doubletfinder --out data/doubletfinder/unified
python scripts/python/doublet_robustness.py --unified data/doubletfinder/unified \
    --out data/doubletfinder/unified --target-cluster 17 \
    --palette data/cluster_colors_18.json --plot-format pdf png
```
Reproduces: overall 4.17% doublets (51,726 → 49,568 cells), cluster 17 OR ≈ 0.36
(depleted, not enriched), cluster-composition r ≈ 0.999.

## Notes
- `unify` flags any cell where `valid > original` (would imply cells created by
  removal) rather than emitting a negative rate; these inputs have none.
- These are the aggregate counts the figure is computed from, not per-cell annotations.
