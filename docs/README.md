# Docs

Supplementary technical documentation that's too detailed for the root
`README.md`/`RUN.md` but doesn't belong in code comments. Two kinds of doc
live here:

1. **A human-readable derivation of a machine-readable source** — read the
   source for automation, read the doc for a reviewer/human.
2. **A worked walkthrough of one script's inputs/usage** — for scripts whose
   CLI has enough moving parts (multiple TSV inputs, a build/genome flag,
   sample-order control) that `--help` alone isn't enough context.

## Current docs

| File | What it is | Source of truth it mirrors |
|------|-----------|------------------------------|
| `methods_parameter_table.md` | Every pinned version/threshold/seed | `../config/params.yaml` — **edit the YAML, not this file, then re-derive** |
| `nexus_cnv_figure.md` | Input schema + full CLI usage for the Nexus WES CNV curve figure | `../scripts/python/plot_nexus_cnv_curves.py` |

## Adding a doc

- If it documents parameters/versions: extend `methods_parameter_table.md` to
  match `params.yaml`'s new key, don't start a second parameter doc.
- If it documents one script's usage: name it `<script_topic>.md`, link it
  from `../FIGURES.md`'s Status column for that panel so it's discoverable
  from the figure map, and keep the CLI example in the doc copy-pasteable
  (i.e. runnable as-is against the checked-in `data/` fixtures where possible).
- Don't duplicate what `../RUN.md` (execution order) or `../FIGURES.md`
  (panel→script map) already own — link to them instead of restating.
