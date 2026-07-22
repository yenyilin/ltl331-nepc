#!/usr/bin/env python3
"""
plot_nexus_cnv_curves.py — curve-style CNV visualisations from the tidy TSVs
emitted by nexus_cnv_loader.py.  Three modes (pick with --mode):

  per_sample (default)
      Small multiples: one step-plot row per sample, sorted by timepoint.
      Within each row the CNV call (-2..+2) is drawn as a filled step curve
      across the genome — red above the y=0 baseline (gain / high gain), blue
      below (loss / big loss).  Makes per-timepoint magnitude visible at a
      glance and complements the integer-coloured heatmap.

  frequency
      Cohort summary on a single panel: at each bin, the fraction of samples
      with a gain plotted upward (red), the fraction with a loss plotted
      downward (blue).  A mirror plot for "how recurrent is this lesion".

  combined
      Stacked: the frequency mirror plot on top, the per-sample small
      multiples below, sharing a single x-axis (chromosome layout).  Best
      when you want both cohort recurrence and per-timepoint detail in one
      figure.

Genome binning + chromosome ordering are reused from plot_nexus_cnv_heatmap.py
so the two figures align on the x-axis 1:1.

Outputs to OUTDIR:
  cnv_curves_per_sample.{pdf,png}     when --mode per_sample
  cnv_curves_frequency.{pdf,png}      when --mode frequency
  cnv_curves_combined.{pdf,png}       when --mode combined
  cnv_curves_matrix.tsv               the n_samples x n_bins integer matrix
                                       used (only when --save-matrix)

Example (per-sample, after running nexus_cnv_loader.py):
  python plot_nexus_cnv_curves.py \\
      --calls-tsv data/nexus_calls.tsv \\
      --descriptors-tsv data/nexus_descriptors.tsv \\
      --bin-size 1000000 --mode per_sample \\
      --outdir nexus_plots --plot-format pdf png
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _utils import setup_style, savefig
from plot_nexus_cnv_heatmap import (
    GENOME_SIZES, CHROM_ORDER, build_bins, paint_sample, sort_samples,
)


GAIN_COLOR = '#B2182B'  # RdBu_r deep red — gain / high gain
LOSS_COLOR = '#2166AC'  # RdBu_r deep blue — loss / big loss
ZERO_LINE  = '#888888'
BUILD = 'hg38'          # genome build label for the x-axis (set in main)


# --------------------------------------------------------------------------- #
def gene_positions(genes_df, bins_df, bin_size):
    """Map each driver gene to its global bin index (by gene-body midpoint).
    Returns [(global_idx, gene_name), ...] ordered along the genome."""
    off = bins_df.groupby('chrom', sort=False)['global_idx'].min().to_dict()
    nb = bins_df.groupby('chrom', sort=False)['global_idx'].size().to_dict()
    pos = []
    for _, r in genes_df.iterrows():
        c = str(r['chrom'])
        if c not in off:
            continue
        mid = (float(r['start']) + float(r['end'])) / 2.0
        b = min(int(mid // bin_size), nb[c] - 1)
        pos.append((off[c] + b, str(r['gene'])))
    return sorted(pos)


def annotate_genes(ax, positions, labels=True):
    """Vertical marker line down the panel; italic gene label above (labels=True
    only on the top panel — pass labels=False to extend the line through the rows)."""
    for gx, name in positions:
        ax.axvline(gx, color='#222222', lw=0.9, ls=(0, (3, 1.4)), alpha=0.9,
                   zorder=6, clip_on=True)
        if labels:
            ax.annotate(name, xy=(gx, 1.0), xycoords=('data', 'axes fraction'),
                        xytext=(0, 11), textcoords='offset points', rotation=90,
                        ha='center', va='bottom', fontsize=6, fontstyle='italic',
                        fontweight='medium', color='#000000', clip_on=False, zorder=7)


# --------------------------------------------------------------------------- #
# zygosity / LOH track (from nexus_zygosity.tsv)
# --------------------------------------------------------------------------- #
ZYG_CODE = {'imbalance': 1, 'hom_high': 2, 'hom_low': 2}   # LOH/homozygous -> 2
ZYG_COLORS = {1: '#762A83', 2: '#B35806'}                  # imbalance=purple, LOH=orange-brown
ZYG_LABELS = {1: 'allelic imbalance', 2: 'LOH'}

# Passages that are the actual scRNA-seq material AND re-measured by WES
# (Supplementary Fig. S3): preCx, 16, 20, 22 wk. Labelled 'Gen11 · …' in the
# sample-order TSV (parallel to Gen3/Gen5/Gen7) and defined in the footnote; the
# '(gen11)' suffix here is only a fallback when no display map is given. The 8/12-wk
# tracks are a separate xenograft generation (LTL331-5) and are NOT flagged.
SCRNA_WES_MATCHED = {'LTL331_preCX1', 'LTL331_16wk1', 'LTL331_20wk', 'LTL331_22wk2'}


def paint_zygosity(segments, bin_size, chrom_offsets, chrom_nbins, n_total):
    """Rasterise one sample's zygosity segments -> length-n_total int vector
    (0 none, 1 imbalance, 2 LOH); LOH takes priority over imbalance in a bin."""
    out = np.zeros(n_total, dtype=np.int8)
    for _, r in segments.iterrows():
        chrom = r['Chromosome']
        if chrom not in chrom_offsets:
            continue
        code = ZYG_CODE.get(str(r['value_label']), 0)
        if code == 0:
            continue
        off, nb = chrom_offsets[chrom], chrom_nbins[chrom]
        b0 = max(0, int(r['Start'] // bin_size))
        b1 = min(nb - 1, int(r['End'] // bin_size))
        seg = out[off + b0: off + b1 + 1]
        np.maximum(seg, code, out=seg)
    return out


def _draw_zyg_band(ax, x, zyg, y0=-3.35, h=0.5):
    """Thin colored band of zygosity states along the bottom of a per-sample row."""
    for code, color in ZYG_COLORS.items():
        m = zyg == code
        if m.any():
            ax.bar(x[m], height=h, bottom=y0, width=1.0, color=color,
                   linewidth=0, align='center', zorder=3)


# --------------------------------------------------------------------------- #
def _chrom_x_axis(bins_df):
    """Return (boundaries [n_chrom+1], centers [n_chrom], labels [n_chrom])
    in global bin-index coordinates."""
    boundaries, centers, labels = [], [], []
    for chrom, idx_series in bins_df.groupby('chrom', sort=False)['global_idx']:
        idxs = idx_series.values
        boundaries.append(idxs[0] - 0.5)
        centers.append(idxs.mean())
        labels.append(chrom.replace('chr', ''))
    boundaries.append(bins_df['global_idx'].max() + 0.5)
    return np.array(boundaries), np.array(centers), labels


def _decorate_chrom_axis(ax, boundaries, centers, labels, show_top_labels,
                        show_xlabel=False):
    for b in boundaries[1:-1]:
        ax.axvline(b, color='#bbbbbb', lw=0.3)
    ax.set_xlim(boundaries[0], boundaries[-1])
    ax.set_xticks([])
    if show_top_labels:
        ax_top = ax.secondary_xaxis('top')
        ax_top.set_xticks(centers)
        ax_top.set_xticklabels(labels, fontsize=8)
        ax_top.tick_params(length=0)
    if show_xlabel:
        ax.set_xlabel(f'chromosome ({BUILD})', fontsize=9)


# --------------------------------------------------------------------------- #
# draw helpers (used by both standalone and combined modes)
# --------------------------------------------------------------------------- #
def _draw_frequency_axes(ax, matrix, x, n_samples, show_legend=True):
    """Mirror gain/loss frequency curve into a pre-made axes."""
    gain_freq = (matrix > 0).mean(axis=0) * 100
    loss_freq = (matrix < 0).mean(axis=0) * 100

    ax.fill_between(x, 0, gain_freq, step='mid',
                    color=GAIN_COLOR, alpha=0.85, linewidth=0,
                    label=f'gain ({n_samples} samples)')
    ax.fill_between(x, 0, -loss_freq, step='mid',
                    color=LOSS_COLOR, alpha=0.85, linewidth=0,
                    label='loss')
    ax.axhline(0, color=ZERO_LINE, lw=0.4)

    ymax = max(50, np.ceil(max(gain_freq.max(), loss_freq.max()) / 10) * 10)
    ax.set_ylim(-ymax, ymax)
    yticks = np.arange(-ymax, ymax + 1, max(10, int(ymax / 5)))
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'{abs(int(t))}%' for t in yticks], fontsize=8)
    ax.set_ylabel('fraction of samples', fontsize=9)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == 'left')
    ax.tick_params(length=0)
    if show_legend:
        ax.legend(loc='upper right', fontsize=8, frameon=False, ncol=2)


def _draw_per_sample_row(ax, row, x, label, zyg=None):
    """Filled step curve for a single sample into a pre-made axes."""
    y = row.astype(float)
    ax.fill_between(x, 0, y, step='mid', where=(y > 0),
                    color=GAIN_COLOR, alpha=0.85, linewidth=0)
    ax.fill_between(x, 0, y, step='mid', where=(y < 0),
                    color=LOSS_COLOR, alpha=0.85, linewidth=0)
    ax.axhline(0, color=ZERO_LINE, lw=0.35)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(['-2', '-1', '0', '+1', '+2'], fontsize=6)
    if zyg is not None:
        ax.set_ylim(-3.7, 2.6)          # headroom for the zygosity band below
        _draw_zyg_band(ax, x, zyg)
    else:
        ax.set_ylim(-2.6, 2.6)
    ax.set_ylabel(label, fontsize=8, rotation=0, ha='right', va='center',
                  labelpad=8)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == 'left')
    ax.tick_params(length=0)


def _sample_label(sample, tp_map, label_map=None):
    if label_map and sample in label_map:
        return label_map[sample]            # display name already encodes the generation (e.g. 'Gen11 · …')
    star = ' (gen11)' if sample in SCRNA_WES_MATCHED else ''   # fallback flag when no display map
    tp = tp_map.get(sample)
    return f'{sample}\n[{tp}]{star}' if isinstance(tp, str) else sample + star


def _tp_map(descr_df):
    if descr_df is None or 'timepoint' not in descr_df.columns:
        return {}
    return descr_df.set_index('sample')['timepoint'].to_dict()


# --------------------------------------------------------------------------- #
# mode 1: per-sample small multiples
# --------------------------------------------------------------------------- #
def plot_per_sample(matrix, sample_order, bins_df, outdir, formats,
                    descr_df=None, title=None, genes=None, zyg=None, label_map=None):
    plt = setup_style(font_size=9)
    n_samples, n_bins = matrix.shape
    x = bins_df['global_idx'].values
    boundaries, centers, labels = _chrom_x_axis(bins_df)
    tp_map = _tp_map(descr_df)

    fig_w = max(10.0, min(20.0, n_bins * 0.004 + 4.0))
    fig_h = max(3.0, n_samples * 0.8 + 1.2)
    fig, axes = plt.subplots(n_samples, 1, figsize=(fig_w, fig_h),
                             sharex=True, sharey=True)
    if n_samples == 1:
        axes = [axes]

    for i, (ax, sample) in enumerate(zip(axes, sample_order)):
        _draw_per_sample_row(ax, matrix[i], x,
                             _sample_label(sample, tp_map, label_map),
                             zyg=(zyg[i] if zyg is not None else None))
        _decorate_chrom_axis(ax, boundaries, centers, labels,
                             show_top_labels=(i == 0),
                             show_xlabel=(i == n_samples - 1))
        if genes:
            annotate_genes(ax, genes, labels=(i == 0))

    # single legend on the top axes
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=GAIN_COLOR, label='gain (+1 / +2)'),
               Patch(facecolor=LOSS_COLOR, label='loss (-1 / -2)')]
    if zyg is not None:
        handles += [Patch(facecolor=ZYG_COLORS[1], label='allelic imbalance'),
                    Patch(facecolor=ZYG_COLORS[2], label='LOH')]
    axes[0].legend(handles=handles, loc='upper right', fontsize=6,
                   frameon=False, ncol=2)
    if title:
        fig.suptitle(title, fontsize=8, y=1.02)

    fig.tight_layout()
    savefig(fig, outdir, 'cnv_curves_per_sample', formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# mode 2: cohort gain/loss frequency mirror plot
# --------------------------------------------------------------------------- #
def plot_frequency(matrix, sample_order, bins_df, outdir, formats,
                   title=None, genes=None):
    plt = setup_style(font_size=9)
    n_samples, n_bins = matrix.shape
    x = bins_df['global_idx'].values
    boundaries, centers, labels = _chrom_x_axis(bins_df)

    fig_w = max(10.0, min(20.0, n_bins * 0.004 + 4.0))
    fig, ax = plt.subplots(figsize=(fig_w, 3.5))
    _draw_frequency_axes(ax, matrix, x, n_samples)
    _decorate_chrom_axis(ax, boundaries, centers, labels,
                         show_top_labels=True, show_xlabel=True)
    if genes:
        annotate_genes(ax, genes)
    if title:
        ax.set_title(title, fontsize=8)
    fig.tight_layout()
    savefig(fig, outdir, 'cnv_curves_frequency', formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# mode 3: combined — frequency top + per-sample below, shared x-axis
# --------------------------------------------------------------------------- #
def plot_combined(matrix, sample_order, bins_df, outdir, formats,
                  descr_df=None, title=None, genes=None, zyg=None, label_map=None,
                  figwidth_mm=None, figheight_mm=None):
    plt = setup_style(font_size=9)
    n_samples, n_bins = matrix.shape
    x = bins_df['global_idx'].values
    boundaries, centers, labels = _chrom_x_axis(bins_df)
    tp_map = _tp_map(descr_df)

    # default authored ~landscape-page size (cap 13 in, not 20) so a width=\textwidth
    # LANDSCAPE placement scales it ~1:1 and the 8-9 pt fonts stay legible.
    # --figwidth-mm / --figheight-mm override exactly (e.g. 255 x 180 for a landscape A4);
    # the 11 rows auto-compress via the gridspec height ratios.
    fig_w = max(10.0, min(13.0, n_bins * 0.004 + 4.0))
    freq_h = 2.5
    per_h_per_row = 0.85   # taller per-sample rows so the -2..+2 CN track breathes
    fig_h = freq_h + n_samples * per_h_per_row + 1.2
    if figwidth_mm:
        fig_w = figwidth_mm / 25.4
    if figheight_mm:
        fig_h = figheight_mm / 25.4

    fig = plt.figure(figsize=(fig_w, fig_h))
    # gridspec: 1 freq panel + n_samples per-sample rows, all same column
    height_ratios = [freq_h] + [per_h_per_row] * n_samples
    gs = fig.add_gridspec(1 + n_samples, 1,
                          height_ratios=height_ratios, hspace=0.18)

    ax_freq = fig.add_subplot(gs[0, 0])
    _draw_frequency_axes(ax_freq, matrix, x, n_samples)
    _decorate_chrom_axis(ax_freq, boundaries, centers, labels,
                         show_top_labels=True, show_xlabel=False)
    # cohort label
    ax_freq.set_title('cohort frequency (top) + \nper-sample CNV (below)',
                      fontsize=9, loc='left', pad=14, x=-0.1)
    if genes:
        annotate_genes(ax_freq, genes, labels=True)

    from matplotlib.patches import Patch
    # per-sample zygosity key: placed on the top panel (upper-left, over the
    # empty high-frequency region) so it never overlaps a per-sample track
    if zyg is not None:
        gain_loss_leg = ax_freq.get_legend()
        zyg_leg = ax_freq.legend(
            handles=[Patch(facecolor=ZYG_COLORS[1], label='allelic imbalance'),
                     Patch(facecolor=ZYG_COLORS[2], label='LOH')],
            loc='upper left', fontsize=8, frameon=False, ncol=2,
            title='per-sample zygosity', title_fontsize=6)
        zyg_leg._legend_box.align = 'left'
        if gain_loss_leg is not None:
            ax_freq.add_artist(gain_loss_leg)
    for i, sample in enumerate(sample_order):
        ax = fig.add_subplot(gs[1 + i, 0], sharex=ax_freq)
        _draw_per_sample_row(ax, matrix[i], x,
                             _sample_label(sample, tp_map, label_map),
                             zyg=(zyg[i] if zyg is not None else None))
        _decorate_chrom_axis(ax, boundaries, centers, labels,
                             show_top_labels=False,
                             show_xlabel=(i == n_samples - 1))
        if genes:
            annotate_genes(ax, genes, labels=False)

    # key for the '*' row-label flag on the scRNA/WES-matched passages
    fig.text(0.008, 0.08, 'Gen11 = the xenograft generation profiled by scRNA-seq and\n'
             're-measured by WES (Fig. S3)',
             fontsize=8, ha='left', va='bottom', style='italic')

    if title:
        fig.suptitle(title, fontsize=10, y=0.995)
    savefig(fig, outdir, 'cnv_curves_combined', formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--calls-tsv', required=True,
                   help='long-format CNV calls (nexus_cnv_loader.py output)')
    p.add_argument('--descriptors-tsv', default=None,
                   help='optional descriptors TSV for timepoint-aware sorting')
    p.add_argument('--bin-size', type=int, default=1_000_000,
                   help='genome bin width in bp (default 1 Mbp; match the '
                        'heatmap script so x-axes align)')
    p.add_argument('--mode',
                   choices=['per_sample', 'frequency', 'combined'],
                   default='per_sample')
    p.add_argument('--outdir', default='nexus_plots')
    p.add_argument('--filter-samples', nargs='+', default=None)
    p.add_argument('--build', choices=['hg38', 'hg19'], default='hg38',
                   help='genome build for chromosome sizes / x-axis (Nexus WGS = hg38)')
    p.add_argument('--genes', default=None,
                   help='TSV of driver loci to annotate on the top panel '
                        '(cols: gene, chrom, start, end; e.g. data/driver_genes_hg38.tsv)')
    p.add_argument('--zygosity', default=None,
                   help='TSV of LOH/allelic-imbalance segments (nexus_zygosity.tsv) — '
                        'drawn as a band below each per-sample CNV row')
    p.add_argument('--sample-order', default=None,
                   help='TSV (cols: sample, [display_name], order) giving the row order '
                        'and optional biological labels; e.g. data/nexus_sample_order.tsv')
    p.add_argument('--title', default=None)
    p.add_argument('--save-matrix', action='store_true',
                   help='also write the painted n_samples x n_bins integer '
                        'matrix as TSV')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    p.add_argument('--figwidth-mm', type=float, default=None,
                   help='combined mode: author the figure at this FINAL width in mm so a '
                        'width=\\textwidth LANDSCAPE placement scales it ~1:1 (e.g. 255 for '
                        'a landscape A4 text column); fonts then stay ~8-9 pt on the page')
    p.add_argument('--figheight-mm', type=float, default=None,
                   help='combined mode: final figure height in mm (e.g. 180 for landscape A4); '
                        'the per-sample rows auto-compress to fit')
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(args.calls_tsv, sep='\t')
    if args.filter_samples:
        calls = calls[calls['sample'].isin(args.filter_samples)]
    if calls.empty:
        raise SystemExit('no calls after filter')

    descr = None
    if args.descriptors_tsv:
        descr = pd.read_csv(args.descriptors_tsv, sep='\t')

    print(f'Loaded {len(calls)} segments across {calls["sample"].nunique()} samples')

    global BUILD
    BUILD = args.build
    bins_df = build_bins(args.bin_size, GENOME_SIZES[args.build])
    n_total = len(bins_df)
    chrom_offsets, chrom_nbins = {}, {}
    for chrom, group in bins_df.groupby('chrom', sort=False):
        chrom_offsets[chrom] = int(group['global_idx'].min())
        chrom_nbins[chrom] = len(group)
    print(f'Built {n_total} bins of {args.bin_size:,} bp across '
          f'{len(chrom_offsets)} chromosomes')

    present = calls['sample'].unique().tolist()
    label_map = None
    if args.sample_order:
        so = pd.read_csv(args.sample_order, sep='\t', comment='#')
        if 'order' in so.columns:
            so = so.sort_values('order')
        wanted = [s for s in so['sample'] if s in present]
        extra = [s for s in present if s not in set(wanted)]
        if extra:
            print(f'[warn] samples not in --sample-order (appended): {extra}')
        sample_order = wanted + sort_samples(extra, descr)
        if 'display_name' in so.columns:
            label_map = dict(zip(so['sample'], so['display_name'].astype(str)))
    else:
        sample_order = sort_samples(present, descr)

    matrix = np.zeros((len(sample_order), n_total), dtype=np.int8)
    for i, s in enumerate(sample_order):
        matrix[i] = paint_sample(calls[calls['sample'] == s], args.bin_size,
                                  chrom_offsets, chrom_nbins, n_total)

    zyg = None
    if args.zygosity:
        zc = pd.read_csv(args.zygosity, sep='\t')
        if args.filter_samples:
            zc = zc[zc['sample'].isin(args.filter_samples)]
        zyg = np.zeros((len(sample_order), n_total), dtype=np.int8)
        for i, s in enumerate(sample_order):
            zyg[i] = paint_zygosity(zc[zc['sample'] == s], args.bin_size,
                                    chrom_offsets, chrom_nbins, n_total)
        print(f'Loaded zygosity for {(zyg != 0).any(axis=1).sum()} samples')

    if args.save_matrix:
        mat_df = pd.DataFrame(matrix, index=sample_order,
                              columns=bins_df['global_idx'].astype(str).values)
        mat_df.to_csv(outdir / 'cnv_curves_matrix.tsv', sep='\t')

    gpos = []
    if args.genes:
        genes_df = pd.read_csv(args.genes, sep='\t', comment='#')
        gpos = gene_positions(genes_df, bins_df, args.bin_size)
        print(f'Annotating {len(gpos)} driver genes on the top panel')

    if args.mode == 'per_sample':
        plot_per_sample(matrix, sample_order, bins_df, outdir, args.plot_format,
                        descr_df=descr, title=args.title, genes=gpos,
                        zyg=zyg, label_map=label_map)
    elif args.mode == 'frequency':
        plot_frequency(matrix, sample_order, bins_df, outdir, args.plot_format,
                       title=args.title, genes=gpos)
    else:  # combined
        plot_combined(matrix, sample_order, bins_df, outdir, args.plot_format,
                      descr_df=descr, title=args.title, genes=gpos,
                      zyg=zyg, label_map=label_map,
                      figwidth_mm=args.figwidth_mm, figheight_mm=args.figheight_mm)
    print('Done.')


if __name__ == '__main__':
    main()
