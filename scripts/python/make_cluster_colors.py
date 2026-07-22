#!/usr/bin/env python3
"""
make_cluster_colors.py — derive an 18-cluster palette where every cluster has a DISTINCT
colour that is a clear SHADE of its phenotype-band colour. So Fig 2A (per-cluster) and the
phenotype-banded panels (plot_timepoint_panel.py --by-phenotype) read as the same colour
language: same hue family = same phenotype, different shade = different cluster.

Method (curated v2, 2026-07-11): each phenotype family is a perceptually-even ramp
interpolated in CIELab between a hand-tuned LIGHT and DARK anchor (see ANCHORS), members
ordered by canonical cluster order (light -> dark). Lab (not HLS/RGB) keeps the steps
visually even, and the anchors are widened so 6-7 same-hue shades stay distinguishable
(min pairwise ΔE ≈ 9) while remaining unmistakably in-family. cl17 (gateway) is a deep
fuchsia singleton — a colourblind-safe pop that also matches the graphical-abstract palette
(data/cluster_colors_gateway.json). Design + colourblind stress-test: see
scratchpad palette18_design*.py history.

Writes <out> (cluster -> hex, in canonical order) and prints the table.

Example:
  python make_cluster_colors.py --labels data/cluster_phenotype_labels.json \\
      --order data/cluster_order.json --out data/cluster_colors_18.json
"""
import argparse, json
from pathlib import Path

import numpy as np

# Curated v2 design (self-contained on purpose: the palette is a fixed, hand-tuned artefact,
# so member->shade order is embedded here rather than read from the drifting cluster_order.json —
# guarantees a re-run reproduces the committed palette byte-for-byte). Each family is a CIELab
# ramp from a LIGHT to a DARK anchor across its members, listed light -> dark.
GROUPS = [
    # (phenotype, [cluster ids: light -> dark], light_anchor, dark_anchor)
    ('AR-high PRAD',       [11, 15],              '#a99bec', '#3f2588'),   # violet
    ('AR-low PRAD',        [0, 2, 5, 8, 12, 16],  '#ff8347', '#6e3603'),   # orange -> deep brown
    ('Pre-EMT',            [4, 6],                '#ffd524', '#b98600'),   # gold
    ('Intermediate (EMT)', [17],                  '#b5179e', '#b5179e'),   # fuchsia singleton
    ('NEPC', [10, 13, 1, 3, 9, 14, 7],            '#5bf58f', '#0e454c'),   # green -> teal (ASCL1+ -> ASCL1-)
]
# emit keys in the canonical file order (matches the historical committed JSON for clean diffs)
OUTPUT_ORDER = [11, 15, 0, 2, 5, 8, 12, 16, 4, 6, 17, 10, 13, 1, 3, 9, 14, 7]


# ---- sRGB <-> CIELab (D65) ----
_M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
_Mi = np.linalg.inv(_M)
_WHITE = np.array([0.95047, 1.0, 1.08883])


def _lin(c):
    c = c / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)


def _hex2rgb(h):
    h = h.lstrip('#'); return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


def _rgb2hex(rgb):
    r = np.clip(np.round(rgb), 0, 255).astype(int); return '#%02x%02x%02x' % (r[0], r[1], r[2])


def _f(t):
    d = 6 / 29; return np.where(t > d ** 3, np.cbrt(t), t / (3 * d * d) + 4 / 29)


def _fi(t):
    d = 6 / 29; return np.where(t > d, t ** 3, 3 * d * d * (t - 4 / 29))


def _rgb2lab(h):
    xyz = _M @ _lin(_hex2rgb(h)); f = _f(xyz / _WHITE)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _lab2hex(lab):
    fy = (lab[0] + 16) / 116; fx = fy + lab[1] / 500; fz = fy - lab[2] / 200
    return _rgb2hex(_srgb(_Mi @ (_WHITE * np.array([_fi(fx), _fi(fy), _fi(fz)]))) * 255.0)


def shades(lo, hi, k):
    """k in-family shades: CIELab ramp from light anchor to dark anchor."""
    if k == 1:
        return [lo]
    a, b = _rgb2lab(lo), _rgb2lab(hi)
    return [_lab2hex(a + (b - a) * i / (k - 1)) for i in range(k)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='data/cluster_colors_18.json')
    args = ap.parse_args()

    colors, pheno_of, anchors_of = {}, {}, {}
    for name, ids, lo, hi in GROUPS:
        for cid, hx in zip(ids, shades(lo, hi, len(ids))):
            colors[cid] = hx; pheno_of[cid] = name; anchors_of[cid] = (lo, hi)

    out = {str(c): colors[c] for c in OUTPUT_ORDER}
    Path(args.out).write_text(json.dumps(out, indent=0) + "\n")

    print(f"{'cl':>3}  {'phenotype':20} {'hex':9}  (anchors light..dark)")
    for c in OUTPUT_ORDER:
        lo, hi = anchors_of[c]
        print(f"{c:>3}  {pheno_of[c]:20} {colors[c]}  ({lo}..{hi})")
    print(f"\nwrote {args.out}  (n={len(out)})")


if __name__ == '__main__':
    main()
