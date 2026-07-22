#!/usr/bin/env python3
"""
plot_nexus_cnv_heatmap.py — sample-by-genome CNV heatmap from the tidy TSVs
emitted by nexus_cnv_loader.py.

Rows: samples, sorted by timepoint (preCx -> 22wk; samples without a parseable
      timepoint go to the bottom).
Cols: fixed-width genome bins of --bin-size bp (default 1 Mbp) tiled across
      chr1..chr22, chrX, chrY (hg19 sizes; manuscript build).
Color: integer CNV call per bin chosen as the max-|value| of overlapping
       segments (so a focal high-gain isn't diluted by neighbouring losses).
         -2 big loss     -1 loss     0 neutral     +1 gain     +2 high gain

Outputs to OUTDIR:
  cnv_heatmap.{pdf,png}      the heatmap (publication-grade, fonttype 42)
  cnv_heatmap_matrix.tsv     the underlying n_samples x n_bins integer matrix
  cnv_heatmap_bins.tsv       global_idx <-> (chrom, start, end) reference

Example:
  python plot_nexus_cnv_heatmap.py \\
      --calls-tsv data/nexus_calls.tsv \\
      --descriptors-tsv data/nexus_descriptors.tsv \\
      --bin-size 1000000 \\
      --outdir nexus_plots --plot-format pdf png
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _utils import setup_style, savefig, timepoint_key


# hg19 chromosome sizes (UCSC). Manuscript build is NCBI Build 37 == hg19.
HG19_SIZES = {
    'chr1':  249250621, 'chr2':  243199373, 'chr3':  198022430,
    'chr4':  191154276, 'chr5':  180915260, 'chr6':  171115067,
    'chr7':  159138663, 'chr8':  146364022, 'chr9':  141213431,
    'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
    'chr16':  90354753, 'chr17':  81195210, 'chr18':  78077248,
    'chr19':  59128983, 'chr20':  63025520, 'chr21':  48129895,
    'chr22':  51304566, 'chrX':  155270560, 'chrY':   59373566,
}
# GRCh38/hg38 primary-assembly chromosome sizes (the Nexus WGS export is Build 38 —
# nexus/927-331.WGS/organism.txt: build=NCBI Build 38). Use these for the Nexus CNV plots.
HG38_SIZES = {
    'chr1':  248956422, 'chr2':  242193529, 'chr3':  198295559,
    'chr4':  190214555, 'chr5':  181538259, 'chr6':  170805979,
    'chr7':  159345973, 'chr8':  145138636, 'chr9':  138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16':  90338345, 'chr17':  83257441, 'chr18':  80373285,
    'chr19':  58617616, 'chr20':  64444167, 'chr21':  46709983,
    'chr22':  50818468, 'chrX':  156040895, 'chrY':   57227415,
}
GENOME_SIZES = {'hg19': HG19_SIZES, 'hg38': HG38_SIZES}
CHROM_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# diverging palette anchored on ColorBrewer RdBu_r, 5 discrete levels
DISCRETE_COLORS = {
    -2: '#2166AC',  # big_loss   (deep blue)
    -1: '#92C5DE',  # loss       (light blue)
     0: '#F7F7F7',  # neutral    (near-white)
     1: '#F4A582',  # gain       (light red)
     2: '#B2182B',  # high_gain  (deep red)
}
LEVELS = [-2, -1, 0, 1, 2]


def _load_label_map(raw):
    """sample id -> y-axis label, from an inline JSON string or a path to one."""
    if raw is None:
        return {}
    import json
    try:
        d = json.loads(raw)
    except json.JSONDecodeError:
        d = json.loads(Path(raw).read_text())
    return {str(k): str(v) for k, v in d.items()}


# --------------------------------------------------------------------------- #
# binning
# --------------------------------------------------------------------------- #
def build_bins(bin_size, chrom_sizes=HG19_SIZES):
    """Tile the genome into right-open bins; returns a long DataFrame."""
    rows, g = [], 0
    for chrom in CHROM_ORDER:
        size = chrom_sizes.get(chrom)
        if size is None:
            continue
        n = int(np.ceil(size / bin_size))
        for i in range(n):
            rows.append({'chrom': chrom, 'bin_idx': i,
                         'start': i * bin_size,
                         'end': min((i + 1) * bin_size, size),
                         'global_idx': g})
            g += 1
    return pd.DataFrame(rows)


def paint_sample(segments, bin_size, chrom_offsets, chrom_nbins, n_total):
    """Rasterise one sample's segments into a length-n_total int vector.

    Resolves overlap by max-|value| (preserves focal high-gain/loss signal).
    """
    out = np.zeros(n_total, dtype=np.int8)
    out_abs = np.zeros(n_total, dtype=np.int8)
    seg = segments[['Chromosome', 'Start', 'End', 'Value']].to_numpy()
    for chrom, s, e, v in seg:
        if chrom not in chrom_offsets:
            continue
        offset = chrom_offsets[chrom]
        n_bins = chrom_nbins[chrom]
        b0 = max(0, int(s) // bin_size)
        b1 = min(n_bins - 1, (int(e) - 1) // bin_size)
        if b1 < b0:
            continue
        mag = abs(int(v))
        slc = slice(offset + b0, offset + b1 + 1)
        mask = mag > out_abs[slc]
        if mask.any():
            sub = out[slc]
            sub_abs = out_abs[slc]
            sub[mask] = int(v)
            sub_abs[mask] = mag
            out[slc] = sub
            out_abs[slc] = sub_abs
    return out


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #
def sort_samples(samples, descr_df):
    """Sort by timepoint: 'parent' (baseline tumor) -> preCx -> 22wk; untyped
    samples appended at the end.
    """
    if descr_df is None or descr_df.empty or 'timepoint' not in descr_df.columns:
        return sorted(samples)
    tp = descr_df.set_index('sample')['timepoint'].to_dict()
    def keyfn(s):
        t = tp.get(s)
        if t is None or (isinstance(t, float) and np.isnan(t)):
            return (1e9, s)  # untyped go last
        if isinstance(t, str) and t.lower() in ('parent', 'baseline', 't0'):
            return (-2.0, s)  # before preCx (which sorts to -1.0)
        return timepoint_key(t) + (s,)
    return sorted(samples, key=keyfn)


def plot_heatmap(matrix, sample_order, bins_df, outdir, formats,
                 descr_df=None, title=None, label_map=None, color_map=None):
    plt = setup_style(font_size=7)
    from matplotlib.colors import ListedColormap, BoundaryNorm

    cmap = ListedColormap([DISCRETE_COLORS[v] for v in LEVELS])
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)

    n_samples, n_bins = matrix.shape
    fig_w = max(8.0, min(20.0, n_bins * 0.005 + 3.0))
    fig_h = max(2.5, n_samples * 0.35 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm,
                   interpolation='nearest')

    # y-axis: sample labels; if timepoint known, annotate "sample (timepoint)"
    tp_map = {}
    if descr_df is not None and 'timepoint' in descr_df.columns:
        tp_map = descr_df.set_index('sample')['timepoint'].to_dict()
    label_map = label_map or {}
    ylabels = []
    for s in sample_order:
        if s in label_map:
            ylabels.append(label_map[s])
            continue
        t = tp_map.get(s)
        ylabels.append(f'{s}  [{t}]' if isinstance(t, str) else s)
    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(ylabels, fontsize=6)
    # color y-tick labels by stage (palette chosen distinct from the blue/red
    # CNV semantics so the label color is not misread as a gain/loss cue)
    if color_map:
        for tick, s in zip(ax.get_yticklabels(), sample_order):
            c = color_map.get(s)
            if c:
                tick.set_color(c)
                tick.set_fontweight('bold')

    # x-axis: hide bin ticks; draw chromosome boundaries + labels on top
    ax.set_xticks([])
    chrom_groups = bins_df.groupby('chrom', sort=False)['global_idx']
    boundaries, centers, labels = [], [], []
    for chrom, idx_series in chrom_groups:
        idxs = idx_series.values
        boundaries.append(idxs[0] - 0.5)
        centers.append(idxs.mean())
        labels.append(chrom.replace('chr', ''))
    boundaries.append(bins_df['global_idx'].max() + 0.5)
    for b in boundaries[1:-1]:
        ax.axvline(b, color='#777', lw=0.3)
    ax_top = ax.secondary_xaxis('top')
    ax_top.set_xticks(centers)
    ax_top.set_xticklabels(labels, fontsize=6)
    ax_top.tick_params(length=0)
    ax.set_xlabel('chromosome (hg19)')
    ax.set_xlim(boundaries[0], boundaries[-1])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # discrete colorbar
    cbar = fig.colorbar(im, ax=ax, ticks=LEVELS, fraction=0.025, pad=0.01,
                        shrink=0.6)
    cbar.ax.set_yticklabels(['big loss', 'loss', 'neutral', 'gain', 'high gain'],
                            fontsize=6)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)

    if title:
        ax.set_title(title, fontsize=8)

    fig.tight_layout()
    savefig(fig, outdir, 'cnv_heatmap', formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--calls-tsv', required=True,
                   help='long-format CNV calls (nexus_cnv_loader.py output)')
    p.add_argument('--descriptors-tsv', default=None,
                   help='optional descriptors TSV for timepoint-aware sorting + '
                        'sample labels')
    p.add_argument('--bin-size', type=int, default=1_000_000,
                   help='genome bin width in bp (default 1 Mbp)')
    p.add_argument('--outdir', default='nexus_plots')
    p.add_argument('--filter-samples', nargs='+', default=None,
                   help='only keep these sample names (default: all in TSV)')
    p.add_argument('--title', default=None)
    p.add_argument('--row-order', nargs='+', default=None,
                   help='explicit top-to-bottom row (sample) order; overrides '
                        'the default timepoint sort')
    p.add_argument('--row-label-map', default=None,
                   help='JSON string or path mapping sample id -> y-axis label')
    p.add_argument('--row-color-map', default=None,
                   help='JSON string or path mapping sample id -> y-axis label '
                        'color (hex); used to color rows by phenotype/stage')
    p.add_argument('--extra-empty-samples', nargs='+', default=None,
                   help='sample ids to include as all-neutral rows even with no '
                        'CNV calls (e.g. the inferCNV reference group)')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(args.calls_tsv, sep='\t')
    # 'sample' round-trips through CSV as int64 when every id is a digit (e.g.
    # cluster ids 0..17); force str so --row-order / --row-label-map / isin match.
    calls['sample'] = calls['sample'].astype(str)
    if args.filter_samples:
        calls = calls[calls['sample'].isin([str(s) for s in args.filter_samples])]
    if calls.empty:
        raise SystemExit('no calls after filter')

    descr = None
    if args.descriptors_tsv:
        descr = pd.read_csv(args.descriptors_tsv, sep='\t')
        descr['sample'] = descr['sample'].astype(str)
        if 'timepoint' in descr.columns:
            descr['timepoint'] = descr['timepoint'].astype(str)

    print(f'Loaded {len(calls)} segments across {calls["sample"].nunique()} samples')

    bins_df = build_bins(args.bin_size)
    n_total = len(bins_df)
    chrom_offsets, chrom_nbins = {}, {}
    for chrom, group in bins_df.groupby('chrom', sort=False):
        chrom_offsets[chrom] = int(group['global_idx'].min())
        chrom_nbins[chrom] = len(group)
    print(f'Built {n_total} bins of {args.bin_size:,} bp across {len(chrom_offsets)} chromosomes')

    label_map = _load_label_map(args.row_label_map)
    color_map = _load_label_map(args.row_color_map)
    universe = set(calls['sample'].unique().tolist())
    if args.extra_empty_samples:
        universe |= {str(s) for s in args.extra_empty_samples}

    if args.row_order:
        row_order = [str(s) for s in args.row_order]
        sample_order = [s for s in row_order if s in universe]
        dropped = [s for s in row_order if s not in universe]
        if dropped:
            print(f'[warn] --row-order ids not present, skipped: {dropped}')
        leftover = sort_samples([s for s in universe if s not in set(row_order)],
                                descr)
        if leftover:
            print(f'[warn] ids not in --row-order, appended at bottom: {leftover}')
        sample_order += leftover
    else:
        sample_order = sort_samples(sorted(universe), descr)

    matrix = np.zeros((len(sample_order), n_total), dtype=np.int8)
    for i, s in enumerate(sample_order):
        sub = calls[calls['sample'] == s]
        matrix[i] = paint_sample(sub, args.bin_size,
                                  chrom_offsets, chrom_nbins, n_total)
        print(f'  painted {s}: {(matrix[i] != 0).sum()} non-neutral bins')

    # persist matrix + bin reference for downstream / re-plotting
    bins_df.to_csv(outdir / 'cnv_heatmap_bins.tsv', sep='\t', index=False)
    mat_df = pd.DataFrame(matrix, index=sample_order,
                          columns=bins_df['global_idx'].astype(str).values)
    mat_df.to_csv(outdir / 'cnv_heatmap_matrix.tsv', sep='\t')

    plot_heatmap(matrix, sample_order, bins_df, outdir, args.plot_format,
                 descr_df=descr, title=args.title, label_map=label_map,
                 color_map=color_map)
    print('Done.')


if __name__ == '__main__':
    main()
