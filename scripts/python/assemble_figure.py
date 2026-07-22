#!/usr/bin/env python3
"""
assemble_figure.py — compose label-free panel PDFs into one multi-panel figure,
vector-preserving, with the layout given as a subplot-mosaic string.

Why a mosaic (not --ncols): 5- and 6-panel figures don't tile cleanly and panels have
mismatched aspect ratios. A mosaic lets a panel span cells, so any count/shape works:

    "ABC                 3x2, six equal panels
     DEF"

    "AAABBB              5 panels: A,B half-width on top; C,D,E third-width below
     CCDDEE"             (6-col grid = LCM of 2 and 3; letters span)

    "AABBCC              5 panels: 3 across top, 2 centered below
     .DDEE."             ('.' = empty cell)

Each distinct letter is one panel (its cells must form a rectangle); the letter is also
its drawn label (source panels carry NO labels — that's the point, for reuse/flexibility).

Squeeze levers:
  --trim                crop each source to its ink bounding box (removes savefig margin) — the
                        biggest space win; measured by a low-dpi render, output stays vector.
  --col-weights / --row-weights   non-uniform tracks for heterogeneous aspects.
  auto row height       if --row-weights/--height omitted, each row's height is derived from its
                        panels' aspect ratios so letterboxing is minimized.

Fit is aspect-preserving (letterbox); --align sets where a panel sits in a larger cell.
--label-offset "dx,dy" nudges every label (negative = outward: left/up, e.g. into the margin).

RECIPE FILE (--layout-file) — a whole figure as one version-controlled file. Lines:
  '#' comment ; '@key: value' directive (any long-option name) ; 'LETTER = path[:page]' panel ;
  anything else = a mosaic row. CLI flags override @directives; --map extends the file map.

Examples:
  python scripts/assemble_figure.py --layout-file scripts/figure_layouts/fig5_example.layout
  python scripts/assemble_figure.py --width double --trim --label-offset=-5,-3 \\
    --layout "AAABBB
              CCDDEE" \\
    --map A=figures/fig5A.pdf B=figures/fig5B.pdf C=figures/fig5C.pdf \\
          D=figures/fig5D.pdf E=figures/fig5E.pdf --out figures/Fig5.pdf --png 600
"""
import argparse

import fitz  # PyMuPDF

MM = 72.0 / 25.4
PRESETS = {"single": 90.0, "1.5col": 140.0, "onehalf": 140.0, "double": 180.0}  # mm

# resolved-setting defaults (CLI value, else @directive, else this)
DEF = dict(width="double", height=None, max_height=230.0, col_weights=None, row_weights=None,
           gap=3.0, margin=2.0, trim=False, trim_dpi=100, trim_thresh=248, align="center",
           no_labels=False, lowercase=False, label_size=11.0, label_font="hebo",
           label_offset="0,0", png=0)
TYPES = dict(height=float, max_height=float, gap=float, margin=float, label_size=float,
             png=int, trim_dpi=int, trim_thresh=int, trim=bool, no_labels=bool, lowercase=bool)


def coerce(key, val):
    t = TYPES.get(key, str)
    if t is bool:
        return str(val).strip().lower() in ("1", "true", "yes", "on")
    return t(val)


def parse_layout(text):
    """mosaic string -> (grid, nrows, ncols, spans{letter:(r0,r1,c0,c1)}); validates rectangles."""
    lines = [ln for ln in text.splitlines() if ln.strip() != ""]
    if not lines:
        raise SystemExit("[abort] empty layout")
    ncols = max(len(ln) for ln in lines)
    grid = [ln.ljust(ncols, ".") for ln in lines]
    cells = {}
    for r, ln in enumerate(grid):
        for c, ch in enumerate(ln):
            if ch != ".":
                cells.setdefault(ch, []).append((r, c))
    spans = {}
    for ch, rcs in cells.items():
        rs = [r for r, _ in rcs]; cs = [c for _, c in rcs]
        r0, r1, c0, c1 = min(rs), max(rs), min(cs), max(cs)
        if set(rcs) != {(r, c) for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)}:
            raise SystemExit(f"[abort] panel '{ch}' is not a rectangle in the mosaic")
        spans[ch] = (r0, r1, c0, c1)
    return grid, len(grid), ncols, spans


def parse_recipe(text):
    """recipe file -> (mosaic_text, map{letter:val}, directives{key:val})."""
    mosaic, mp, directives = [], {}, {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        s = s.split("#", 1)[0].strip()        # strip inline comment ('#' not used in paths)
        if not s:
            continue
        if s.startswith("@"):
            body = s[1:]
            sep = ":" if ":" in body else ("=" if "=" in body else None)
            if not sep:
                raise SystemExit(f"[abort] bad directive (need '@key: value'): {line}")
            k, v = body.split(sep, 1)
            directives[k.strip().replace("-", "_")] = v.strip()
        elif "=" in s:
            letter, val = s.split("=", 1)
            mp[letter.strip()] = val.strip()
        else:
            mosaic.append(s)            # '.' preserved (not whitespace)
    return "\n".join(mosaic), mp, directives


def parse_map_pairs(pairs, into):
    for p in pairs or []:
        if "=" not in p:
            raise SystemExit(f"[abort] --map entry '{p}' must be LETTER=path[:page]")
        letter, val = p.split("=", 1)
        into[letter.strip()] = val.strip()
    return into


def resolve_path_page(val):
    page = 0
    if ":" in val and val.rsplit(":", 1)[1].isdigit():
        val, pg = val.rsplit(":", 1); page = int(pg) - 1
    return val, page


def weights(arg, n):
    if not arg:
        return [1.0] * n
    w = [float(x) for x in str(arg).split(",")]
    if len(w) != n:
        raise SystemExit(f"[abort] expected {n} weights, got {len(w)}: {arg}")
    return w


def content_clip(page, dpi, thresh):
    """source ink bbox (points) via a low-dpi gray render; output stays vector."""
    import numpy as np
    pix = page.get_pixmap(dpi=dpi, colorspace=fitz.csGRAY)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    mask = a < thresh
    if not mask.any():
        return page.rect
    ys, xs = np.where(mask)
    s = 72.0 / dpi
    r = fitz.Rect(xs.min() * s, ys.min() * s, (xs.max() + 1) * s, (ys.max() + 1) * s)
    r += (page.rect.x0, page.rect.y0, page.rect.x0, page.rect.y0)
    return r & page.rect


def panel_min_fontsize(page, clip=None):
    """Smallest vector-text font size (pt, source coords) on the page, optionally
    restricted to text intersecting `clip`. Returns None if the page has no
    extractable text (outlined/rasterized panels)."""
    sizes = []
    for b in page.get_text("dict").get("blocks", []):
        for ln in b.get("lines", []):
            for sp in ln.get("spans", []):
                if not sp.get("text", "").strip():
                    continue
                if clip is not None and not fitz.Rect(sp["bbox"]).intersects(clip):
                    continue
                if sp.get("size", 0) > 0:
                    sizes.append(sp["size"])
    return min(sizes) if sizes else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layout", help="mosaic string")
    ap.add_argument("--layout-file", help="recipe file: mosaic + 'LETTER=path' + '@key: value'")
    ap.add_argument("--map", nargs="+", help="LETTER=path[:page] ... (extends a recipe's map)")
    ap.add_argument("--out", help="output PDF (or @out: in the recipe)")
    # overridable settings: default None so 'explicitly set on CLI' is detectable
    ap.add_argument("--width", default=None, help="preset (single/1.5col/double) or mm number")
    ap.add_argument("--height", type=float, default=None, help="canvas height mm (else auto)")
    ap.add_argument("--max-height", type=float, default=None, help="cap on auto height (mm)")
    ap.add_argument("--col-weights", default=None)
    ap.add_argument("--row-weights", default=None)
    ap.add_argument("--gap", type=float, default=None, help="inter-panel gap (mm)")
    ap.add_argument("--margin", type=float, default=None, help="outer margin (mm)")
    ap.add_argument("--trim", action="store_const", const=True, default=None,
                    help="crop each panel to its ink bbox")
    ap.add_argument("--trim-dpi", type=int, default=None)
    ap.add_argument("--trim-thresh", type=int, default=None)
    ap.add_argument("--align", default=None, choices=[None, "center", "topleft", "top", "left"])
    ap.add_argument("--no-labels", action="store_const", const=True, default=None)
    ap.add_argument("--lowercase", action="store_const", const=True, default=None)
    ap.add_argument("--label-size", type=float, default=None)
    ap.add_argument("--label-font", default=None, help="fitz base-14 font (hebo=Helvetica-Bold)")
    ap.add_argument("--label-offset", default=None,
                    help='"dx,dy" in pt; negative = outward (left/up), e.g. "-5,-3" for margin labels')
    ap.add_argument("--png", type=int, default=None, help="also write a PNG at this DPI")
    ap.add_argument("--qa", action="store_true",
                    help="print per-panel scale + effective on-page min font size")
    ap.add_argument("--min-font", type=float, default=6.0,
                    help="QA threshold (pt): warn when a panel's smallest text falls below this on the page")
    args = ap.parse_args()

    # ---- gather recipe (file) then resolve CLI > @directive > DEF ----
    rec_mosaic, rec_map, rec_dir = "", {}, {}
    if args.layout_file:
        rec_mosaic, rec_map, rec_dir = parse_recipe(open(args.layout_file).read())

    def setting(key):
        v = getattr(args, key)
        if v is not None:
            return v
        if key in rec_dir:
            return coerce(key, rec_dir[key])
        return DEF[key]

    text = args.layout or rec_mosaic
    if not text:
        raise SystemExit("[abort] no mosaic (use --layout or a recipe with mosaic rows)")
    out_path = args.out or rec_dir.get("out")
    if not out_path:
        raise SystemExit("[abort] no output (use --out or '@out: …' in the recipe)")

    panel_src = dict(rec_map)
    parse_map_pairs(args.map, panel_src)               # CLI --map overrides per letter
    if not panel_src:
        raise SystemExit("[abort] no panel map (recipe 'LETTER=path' lines or --map)")

    grid, nrows, ncols, spans = parse_layout(text)
    missing = [ch for ch in spans if ch not in panel_src]
    if missing:
        raise SystemExit(f"[abort] no path for panel(s): {missing}")

    width = setting("width")
    W = (PRESETS[width] if width in PRESETS else float(width)) * MM
    gap, margin = setting("gap") * MM, setting("margin") * MM
    trim = setting("trim")
    dx_lab, dy_lab = (float(x) for x in str(setting("label_offset")).split(","))

    # column geometry
    cw = weights(setting("col_weights"), ncols)
    avail_w = W - 2 * margin - (ncols - 1) * gap
    col_w = [avail_w * x / sum(cw) for x in cw]
    col_x = [margin + sum(col_w[:j]) + j * gap for j in range(ncols)]

    def span_w_x(c0, c1):
        return sum(col_w[c0:c1 + 1]) + (c1 - c0) * gap, col_x[c0]

    # open sources + clips (clip aspect drives auto row heights)
    src, clip = {}, {}
    for ch in spans:
        path, pg = resolve_path_page(panel_src[ch])
        d = fitz.open(path)
        if pg >= d.page_count:
            raise SystemExit(f"[abort] panel '{ch}': page {pg+1} not in {path}")
        src[ch] = (d, pg)
        clip[ch] = content_clip(d[pg], setting("trim_dpi"), setting("trim_thresh")) if trim else d[pg].rect

    # row geometry
    if setting("row_weights") or setting("height") is not None:
        rw = weights(setting("row_weights"), nrows)
        H = (setting("height") if setting("height") else setting("max_height")) * MM
        avail_h = H - 2 * margin - (nrows - 1) * gap
        row_h = [avail_h * x / sum(rw) for x in rw]
    else:
        row_h = [0.0] * nrows
        for ch, (r0, r1, c0, c1) in spans.items():
            wpt, _ = span_w_x(c0, c1)
            asp = clip[ch].width / clip[ch].height
            nat = (wpt / asp) / (r1 - r0 + 1)
            for r in range(r0, r1 + 1):
                row_h[r] = max(row_h[r], nat)
        H = sum(row_h) + 2 * margin + (nrows - 1) * gap
        cap = setting("max_height") * MM
        if H > cap:
            s = (cap - 2 * margin - (nrows - 1) * gap) / sum(row_h)
            row_h = [h * s for h in row_h]; H = cap
    row_y = [margin + sum(row_h[:i]) + i * gap for i in range(nrows)]

    # render
    out = fitz.open()
    page = out.new_page(width=W, height=H)
    align = setting("align")
    ax = {"center": 0.5, "left": 0.0, "topleft": 0.0, "top": 0.5}[align]
    ay = {"center": 0.5, "left": 0.5, "topleft": 0.0, "top": 0.0}[align]
    no_labels, lowercase = setting("no_labels"), setting("lowercase")
    lsize, lfont = setting("label_size"), setting("label_font")
    qa_rows = []
    for ch in sorted(spans):
        r0, r1, c0, c1 = spans[ch]
        wpt, x = span_w_x(c0, c1)
        hpt = sum(row_h[r0:r1 + 1]) + (r1 - r0) * gap
        cell = fitz.Rect(x, row_y[r0], x + wpt, row_y[r0] + hpt)
        cl = clip[ch]; asp = cl.width / cl.height
        if wpt / hpt > asp:
            h = hpt; w = h * asp
        else:
            w = wpt; h = w / asp
        tx = cell.x0 + ax * (wpt - w); ty = cell.y0 + ay * (hpt - h)
        d, pg = src[ch]
        page.show_pdf_page(fitz.Rect(tx, ty, tx + w, ty + h), d, pg, clip=cl)
        qa_rows.append((ch, cl.width, w, panel_min_fontsize(d[pg], cl)))
        if not no_labels:
            lab = ch.lower() if lowercase else ch.upper()
            page.insert_text((cell.x0 + dx_lab, cell.y0 + lsize + dy_lab), lab,
                             fontsize=lsize, fontname=lfont, color=(0, 0, 0))

    out.save(out_path, deflate=True)
    print(f"[ok] {out_path}  ({W/MM:.0f}x{H/MM:.0f} mm, {len(spans)} panels, "
          f"trim={'on' if trim else 'off'}, label-offset {dx_lab:g},{dy_lab:g})")

    # ---- QA: effective on-page font sizes (scale = rendered_w / source_w) ----
    thr, flagged = args.min_font, []
    if args.qa:
        print(f"  QA  panel  scale  onpage_mm  src_min_pt  eff_min_pt  (threshold {thr:g}pt; "
              f"author each panel at onpage_mm width so scale->1)")
    for ch, srcw, w, mf in sorted(qa_rows):
        sc = (w / srcw) if srcw else float("nan")
        eff = (mf * sc) if mf else None
        if eff is not None and eff < thr:
            flagged.append((ch, eff))
        if args.qa:
            print(f"  QA  {ch:<5}  {sc:5.2f}  {w / MM:8.1f}  "
                  f"{('%.1f' % mf) if mf else 'n/a':>9}  "
                  f"{('%.1f' % eff) if eff is not None else 'n/a':>9}"
                  f"{'   <-- below threshold' if (eff is not None and eff < thr) else ''}")
    if flagged:
        worst = ", ".join(f"{c}={e:.1f}pt" for c, e in sorted(flagged, key=lambda x: x[1]))
        print(f"[warn] {len(flagged)} panel(s) below {thr:g}pt on the page after scaling: {worst}")
    png = setting("png")
    if png:
        page.get_pixmap(dpi=png).save(out_path.rsplit(".", 1)[0] + ".png")
        print(f"[ok] {out_path.rsplit('.',1)[0]}.png ({png} dpi)")


if __name__ == "__main__":
    main()
