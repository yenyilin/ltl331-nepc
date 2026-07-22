#!/usr/bin/env python3
"""
plot_enrichment_violin.py — Fig 2D: per-cluster paired PRAD vs NE signature-score violins,
clusters on x in the canonical adeno->NE order so the two scores CROSS OVER left->right
(the transdifferentiation), and the non-significant clusters land in the middle (the
intermediate) — which directly answers "how do 2D clusters reflect the longitudinal
transdifferentiation?".

Per cluster: blue = PRAD signature, red = NE signature, with a paired Wilcoxon test
(PRAD vs NE across the cluster's cells) shown as stars (**** <1e-4, *** <1e-3, ** <1e-2,
* <0.05, ns otherwise).

Score source (per side, pick one):
  --prad-key / --ne-key      obs column with a precomputed score (default UCell columns)
  --prad-genes / --ne-genes  a gene set (file: one gene per line, or comma-separated) ->
                             scored with sc.tl.score_genes(use_raw=True) (the .raw lognorm,
                             since .X is scaled). Use this to define "PRAD signature"
                             explicitly (addresses "are PRAD markers = AR?").

Large-n caveat: with thousands of cells per cluster a paired test is significant for almost
any difference. --min-effect lets you require |median(PRAD)-median(NE)| >= that value in
ADDITION to p<0.05 before calling a cluster significant, so "ns" robustly marks clusters
where the two programs are genuinely indistinguishable (the intermediate). Default 0.

Example:
  python plot_enrichment_violin.py --h5ad ltl331_base.h5ad \\
      --cluster-order data/cluster_order.json --cluster-colors data/cluster_colors_18.json \\
      --prad-key AR_pathway_UCell --ne-key Neuroendocrine_UCell \\
      --min-effect 0.03 --outdir figures --stem fig2D
"""
import argparse, json
from pathlib import Path

import numpy as np


def stars(p):
    if p is None or not np.isfinite(p): return 'na'
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'ns'


def load_gene_list(spec):
    p = Path(spec)
    if p.exists():
        toks = []
        for line in p.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            # accept "<gene>" or "<group>\t<gene>" or comma lists
            parts = line.replace(',', '\t').split('\t')
            toks.append(parts[-1].strip())
        return [t for t in toks if t]
    return [g.strip() for g in spec.replace(',', ' ').split() if g.strip()]


def get_score(adata, key, genes, name):
    """Return (per-cell score array, label). Prefer obs key; else score a gene set from raw."""
    import scanpy as sc
    if key and key in adata.obs:
        return np.asarray(adata.obs[key].values, float), key
    if genes:
        gl = load_gene_list(genes)
        src_var = set(map(str, (adata.raw.var_names if adata.raw is not None else adata.var_names)))
        present = [g for g in gl if g in src_var]
        missing = [g for g in gl if g not in src_var]
        if missing:
            print(f'[warn] {name}: {len(missing)} genes absent and skipped: {missing}')
        if not present:
            raise SystemExit(f'{name}: no genes resolve in the object')
        sc.tl.score_genes(adata, present, score_name=name, use_raw=(adata.raw is not None))
        return np.asarray(adata.obs[name].values, float), name
    raise SystemExit(f'{name}: provide an obs key or a gene set')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--cluster-order', default='data/cluster_order.json')
    ap.add_argument('--cluster-colors', default=None)
    ap.add_argument('--prad-key', default='AR_pathway_UCell')
    ap.add_argument('--ne-key', default='Neuroendocrine_UCell')
    ap.add_argument('--prad-genes', default=None)
    ap.add_argument('--ne-genes', default=None)
    ap.add_argument('--prad-label', default='PRAD/AR-pathway score')
    ap.add_argument('--ne-label', default='NE signature')
    ap.add_argument('--min-effect', type=float, default=0.0,
                    help='require |median(PRAD)-median(NE)| >= this AND p<0.05 to call significant')
    ap.add_argument('--exclude-clusters', nargs='*', default=['18'])
    ap.add_argument('--prad-color', default='#2c6fbb')
    ap.add_argument('--ne-color', default='#c0392b')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--stem', default='fig2D')
    ap.add_argument('--figwidth-mm', type=float, default=None,
                    help='author the panel at this FINAL width in mm (the assemble_figure --qa '
                         'onpage_mm target) so the assembler scales it ~1:1 and fonts stay readable; '
                         'note the saved PDF is wider than this by the outside legend — tune with --qa')
    ap.add_argument('--figheight-mm', type=float, default=None, help='panel height in mm (default 4.6 in)')
    ap.add_argument('--no-stars', action='store_true',
                    help='hide the per-cluster significance stars AND the footnote; put the '
                         'paired-Wilcoxon result in the figure caption instead (declutters 2E, '
                         'frees vertical room so the violins fill the panel)')
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from scipy.stats import wilcoxon
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    if args.exclude_clusters:
        keep = ~adata.obs[args.cluster_key].astype(str).isin([str(c) for c in args.exclude_clusters])
        if (~keep).any():
            adata = adata[keep.values].copy()
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)

    prad, prad_src = get_score(adata, args.prad_key, args.prad_genes, 'PRAD_score')
    ne, ne_src = get_score(adata, args.ne_key, args.ne_genes, 'NE_score')
    print(f'PRAD score from {prad_src!r}; NE score from {ne_src!r}')

    cl = adata.obs[args.cluster_key].values
    cats = sorted(set(cl), key=lambda s: int(s) if s.isdigit() else 1e9)
    if Path(args.cluster_order).exists():
        co = [str(c) for c in json.load(open(args.cluster_order)).get('order', [])]
        order = [c for c in co if c in cats] + [c for c in cats if c not in co]
    else:
        order = cats
    colors = {str(k): v for k, v in json.load(open(args.cluster_colors)).items()} \
        if (args.cluster_colors and Path(args.cluster_colors).exists()) else {}

    # per-cluster data + paired test
    data_p, data_n, sig, tops, statrows = [], [], [], [], []
    for c in order:
        m = cl == c
        p_c, n_c = prad[m], ne[m]
        data_p.append(p_c); data_n.append(n_c)
        pv, eff = float('nan'), float('nan')
        try:
            d = p_c - n_c
            if np.allclose(d, 0):
                star = 'ns'; eff = 0.0
            else:
                _, pv = wilcoxon(p_c, n_c)
                eff = abs(np.median(p_c) - np.median(n_c))
                star = stars(pv) if eff >= args.min_effect else 'ns'
        except ValueError:
            star = 'na'
        sig.append(star)
        statrows.append((c, int(m.sum()),
                         float(np.median(p_c)) if len(p_c) else float('nan'),
                         float(np.median(n_c)) if len(n_c) else float('nan'),
                         eff, pv, star))
        tops.append(max(np.nanmax(p_c) if len(p_c) else 0, np.nanmax(n_c) if len(n_c) else 0))

    # per-cluster stats table so 'ns' can be attributed to p-value vs effect-size gate
    print(f'[stats] min_effect gate = {args.min_effect}')
    print(f'{"cluster":>7} {"n":>6} {"med_PRAD":>9} {"med_NE":>9} {"|dmed|":>8} {"p-value":>11}  star')
    for c, n, mp, mn, e, pv, star in statrows:
        gate = '' if (np.isnan(e) or e >= args.min_effect) else '  <fails effect gate>'
        print(f'{c:>7} {n:>6} {mp:>9.4f} {mn:>9.4f} {e:>8.4f} {pv:>11.3g}  {star}{gate}')

    ncl = len(order)
    _dw, _dh = max(7, 0.6 * ncl + 2), 4.6          # authored default (inches)
    if args.figwidth_mm and args.figheight_mm:
        fw, fh = args.figwidth_mm / 25.4, args.figheight_mm / 25.4
    elif args.figwidth_mm:                          # width only -> preserve authored aspect
        fw = args.figwidth_mm / 25.4; fh = _dh * fw / _dw
    elif args.figheight_mm:
        fh = args.figheight_mm / 25.4; fw = _dw * fh / _dh
    else:
        fw, fh = _dw, _dh
    fig, ax = plt.subplots(figsize=(fw, fh))
    off, w = 0.20, 0.34

    def draw(vals_list, xs, color):
        parts = ax.violinplot(vals_list, positions=xs, widths=w, showmeans=False,
                              showextrema=False, showmedians=True)
        for b in parts['bodies']:
            b.set_facecolor(color); b.set_edgecolor('black'); b.set_linewidth(0.3); b.set_alpha(0.85)
        parts['cmedians'].set_color('black'); parts['cmedians'].set_linewidth(1.0)

    draw(data_p, [i - off for i in range(ncl)], args.prad_color)
    draw(data_n, [i + off for i in range(ncl)], args.ne_color)

    ymax = max(tops); ymin = min(min(np.nanmin(p) for p in data_p),
                                 min(np.nanmin(n) for n in data_n))
    pad = (ymax - ymin) * 0.10
    # stars rotated 90 deg so wide '****' don't collide with adjacent clusters at narrow widths
    if not args.no_stars:
        for i, st in enumerate(sig):
            ax.text(i, ymax + pad * 0.5, st, ha='center', va='bottom', rotation=90,
                    fontsize=7 if st != '****' else 6, fontweight='bold')
    top_pad = 1.3 if args.no_stars else 3.0     # no stars -> less headroom, violins fill more
    ax.set_ylim(ymin - pad * 0.5, ymax + pad * top_pad)
    ax.set_xlim(-0.6, ncl - 0.4)
    ax.set_xticks(range(ncl)); ax.set_xticklabels(order, fontsize=7)
    if colors:
        for t, c in zip(ax.get_xticklabels(), order):
            t.set_color(colors.get(c, 'black'))
    ax.set_xlabel('cluster (adenocarcinoma → neuroendocrine)', fontsize=8)
    ax.set_ylabel('signature enrichment score', fontsize=8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    # lower-center anchor -> legend sits ABOVE the axes and grows upward; lifted to 1.10 so it
    # clears the top violins / stars.
    ax.legend(handles=[Patch(facecolor=args.prad_color, label=args.prad_label),
                       Patch(facecolor=args.ne_color, label=args.ne_label)],
              loc='lower center', bbox_to_anchor=(0.5, 1.10), ncol=2, frameon=False,
              fontsize=6.5, handlelength=1.1, handleheight=1.1, columnspacing=1.0,
              handletextpad=0.4, borderpad=0.2)
    if not args.no_stars:
        note = '**** p<1e-4, *** <1e-3, ** <1e-2, * <0.05, ns = not significant (paired Wilcoxon'
        note += f', |Δmedian|≥{args.min_effect})' if args.min_effect else ')'
        ax.text(0.0, -0.55, note, transform=ax.transAxes, fontsize=4, color='0.35')

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.plot_format:
        fig.savefig(outdir / f'{args.stem}.{fmt}', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'[ok] wrote {args.stem}.{{{",".join(args.plot_format)}}}  '
          f'(sig: {dict(zip(order, sig))})')


if __name__ == '__main__':
    main()
