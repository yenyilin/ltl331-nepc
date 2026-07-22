# Assembling multi-panel figures (Fig 2–6) with `assemble_figure.py`

`scripts/assemble_figure.py` is the **final compositing step** for the manuscript
figures. It does not run any analysis — it tiles pre-rendered panel PDFs into one
multi-panel figure, vector-preserving, and draws the panel labels itself.

## The core idea

Each per-panel script (`plot_*.py`) renders **one label-free single-panel PDF**
(no "A"/"B" letters baked in — that is the whole point). `assemble_figure.py` then
mosaics those panels into the final figure and draws the A/B/C… labels. This keeps
panels reusable and lets you re-letter or re-order without re-running the analysis.

## Two ways to run it

### Recommended — a recipe file (one per figure)

Copy `fig5_example.layout`, edit it, and run:

```bash
python scripts/assemble_figure.py --layout-file scripts/figure_layouts/fig5.layout
```

The recipe holds everything — layout, settings, panel map, output path. Line types:
`#` comment | `@key: value` directive | `LETTER = path[:page]` panel | anything else = a mosaic row.

```
@width: double          # single | 1.5col | double, or a number in mm
@trim: true             # crop each panel to its ink bbox (kills savefig margins) — biggest space win
@row-weights: 1.25, 1   # top row taller than bottom
@gap: 3
@label-size: 11
@label-offset: -6, -2   # nudge labels up/left into the margin (negative = outward)
@out: figures/Fig5.pdf
@png: 600               # also emit a raster proof at this DPI

AAABBB                  # the mosaic: A,B half-width on top…
CCDDEE                  #             C,D,E third-width below

A = figures/fig5/fig5A_scenic.pdf
B = figures/fig5/fig5B_decoupler.pdf
C = figures/fig5/fig5C_branch.pdf
D = figures/fig5/fig5D_chromatin.pdf
E = figures/fig5/fig5E_extra.pdf
```

### Or inline via CLI (one-off)

Same thing, no file:

```bash
python scripts/assemble_figure.py --width double --trim --label-offset=-5,-3 \
  --layout "AAABBB
            CCDDEE" \
  --map A=figures/fig5A.pdf B=figures/fig5B.pdf C=figures/fig5C.pdf \
        D=figures/fig5D.pdf E=figures/fig5E.pdf \
  --out figures/Fig5.pdf --png 600
```

CLI flags override `@directives`; `--map` extends/overrides the recipe's panel map.
Use `path:page` (1-indexed) to pull a specific page out of a multi-page PDF.

## The mosaic language

Each **distinct letter = one panel**, and that letter is also its drawn label. A
letter's cells must form a rectangle, so panels can span cells — which is how
mismatched panel counts/aspects tile cleanly:

- `ABC / DEF` → 3×2, six equal panels
- `AAABBB / CCDDEE` → 5 panels (A,B half-width top; C,D,E third-width bottom); the 6-col grid is the LCM of 2 and 3
- `AABBCC / .DDEE.` → `.` = empty cell, so the bottom row's 2 panels center under the top row's 3

So for each of Fig 2–6: write one `.layout` with the mosaic that fits that figure's
panels, point each letter at its label-free `plot_*.py` output PDF, and run `--layout-file`.

## Key knobs

- `--trim` — crops each source to its ink bounding box (removes matplotlib's savefig
  margin). Biggest space win; measured by a throwaway low-dpi render, but output stays vector.
- Auto row height — if you omit `--row-weights`/`--height`, each row's height is derived
  from its panels' aspect ratios to minimize letterboxing (capped by `--max-height`, default 230 mm).
- `--width` presets — `single`=90 mm, `1.5col`/`onehalf`=140 mm, `double`=180 mm (journal column widths).
- `--label-offset "dx,dy"` — negative = outward (into the margin); the usual choice so
  labels don't overlap panel content.
- `--png <dpi>` — writes a raster proof alongside the PDF (e.g. 600) for quick visual checks.
- `--no-labels` / `--lowercase` — suppress or lowercase the drawn labels.
- `--align` — where a panel sits inside a larger cell (center / topleft / top / left).

## Getting started

Copy **`template.layout`** (a fully annotated starter with every `@directive` explained),
edit its mosaic and panel map for your figure, and run it:

```bash
python scripts/assemble_figure.py --layout-file your_figure.layout
```

Write one `.layout` per figure. Point each mosaic letter at the label-free PDF that the
corresponding `plot_*.py` script renders, and `assemble_figure.py` composites them and draws
the panel labels.
