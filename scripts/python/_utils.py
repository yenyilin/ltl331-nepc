"""
_utils.py — shared helpers for the LTL331 revision scripts.

setup_style() / savefig() handle publication-grade plot output (fonttype 42,
Arial, Agg backend). numkey() / timepoint_key() are the sort keys used across
multiple scripts. to_dense_1d() papers over sparse vs dense in scanpy slices.

Sibling import: every script in this dir is launched as `python scripts/foo.py`,
which puts `scripts/` on sys.path[0] -- `from _utils import ...` then resolves
without an __init__.py / package setup.
"""
from pathlib import Path

import numpy as np


def mm_to_in(*vals):
    """Convert millimetres to inches (matplotlib figsize is in inches)."""
    out = tuple(v / 25.4 for v in vals)
    return out[0] if len(out) == 1 else out


def figsize_mm(width_mm, height_mm):
    """figsize (inches) for a panel authored at its FINAL published size in mm.

    READABILITY CONTRACT: author each panel at the width of the layout cell it
    will occupy (see scripts/figure_layouts/*.layout), so assemble_figure.py
    scales it ~1:1 and the point sizes below survive to the printed page. A panel
    authored much wider than its cell gets shrunk by the assembler and its text
    becomes unreadable — run `assemble_figure.py --qa` to catch that.
    """
    return (width_mm / 25.4, height_mm / 25.4)


def setup_style(font_size=8, tick_size=None, label_size=None):
    """Configure matplotlib for editable-text vector output at a fixed point size.

    fonttype 42 embeds TrueType so PDF/PS text stays selectable in
    Illustrator/Inkscape (the default Type 3 outlines text and is rejected
    by some journals). Imports matplotlib lazily so a metrics-only path
    has no plotting dep.

    Font sizes are pinned explicitly (not left as matplotlib's relative
    'medium'/'large') so every text element is a known point size at the panel's
    authored scale. Pair with `figsize_mm(...)` and author at the final cell size
    so these points reach the page unscaled. Journal minimum is ~5-7 pt.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tick_size = tick_size if tick_size is not None else font_size
    label_size = label_size if label_size is not None else font_size
    matplotlib.rcParams.update({
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': font_size, 'axes.linewidth': 0.6,
        # pin the full hierarchy so titles/labels/ticks/legend don't drift to
        # matplotlib's larger relative defaults ('large' titles etc.)
        'axes.titlesize': label_size, 'axes.labelsize': label_size,
        'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size,
        'legend.fontsize': tick_size, 'legend.title_fontsize': label_size,
        'figure.titlesize': label_size,
        'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
        # imshow default 'antialiased' renders thin/1-row colour strips BLANK in the
        # PDF vector backend (fine in PNG/Agg). 'nearest' avoids that path everywhere,
        # incl. artists we can't reach (e.g. CellRank heatmap's fate-probability band).
        'image.interpolation': 'nearest',
    })
    return plt


def savefig(fig, outdir, stem, formats):
    """Save one figure to one or more formats; print each path."""
    for fmt in formats:
        p = Path(outdir) / f'{stem}.{fmt}'
        fig.savefig(p, format=fmt)
        print(f'  wrote {p}')


def numkey(c):
    """Sort key: numeric-looking strings first, sorted numerically; rest last."""
    s = str(c)
    return int(s) if s.isdigit() else 10**9


def timepoint_key(s):
    """Sort key for timepoints: 'preCx'/'pre*' first (-1), then by leading digits."""
    s = str(s)
    if 'pre' in s.lower():
        return (-1.0, s)
    digits = ''.join(ch for ch in s if (ch.isdigit() or ch == '.'))
    return (float(digits) if digits else 1e9, s)


def to_dense_1d(x):
    """Coerce a (possibly sparse) 1-row/column array to a 1D dense numpy array."""
    import scipy.sparse as sp
    return x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
