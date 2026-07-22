"""Proliferation enrichment violin per cluster (original Fig 3F).

sc.tl.score_genes (scanpy equivalent of Seurat's AddModuleScore) on a curated
proliferation gene set, then per-cluster violins. Default set = Tirosh G2/M + S
core (matches the cc.genes Seurat uses for cell-cycle scoring). Override with
--gene-list to drop in the exact set used in the original methods.

Outputs (in --outdir):
  - fig3f_proliferation_violin.{pdf,png}
  - fig3f_proliferation_score_summary.tsv  (per-cluster count/mean/median/q25/q75)
  - fig3f_proliferation_genes_used.tsv     (which genes scored, which missing)
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import scanpy as sc

from _utils import setup_style, savefig, numkey


# Tirosh G2/M + S core. Same superset as Seurat::cc.genes (Tirosh 2016).
DEFAULT_PROLIF = [
    # G2/M
    "MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
    "BUB1", "BUB1B", "CENPE", "CENPF", "BIRC5", "CDC20", "PLK1",
    "UBE2C", "HMMR", "NUSAP1", "TPX2", "KIF11", "KIF20A", "KIF23",
    # S phase / replication
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "PCNA", "RRM1", "RRM2",
    "TYMS", "FEN1", "GINS2", "GMNN", "CDC6", "CDC45", "DTL", "ATAD2",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", default="scripts/fig3f")
    p.add_argument("--cluster-col", default="seurat_clusters")
    p.add_argument("--gene-list", default=None,
                   help="Optional file (one gene per line, '#' = comment) overriding "
                        "the default proliferation set.")
    p.add_argument("--score-name", default="proliferation_score")
    p.add_argument("--score-key", default=None,
                   help="Use a PRECOMPUTED per-cell score column already in obs "
                        "(e.g. 'Proliferation_UCell') instead of re-scoring a gene set. "
                        "When set, --gene-list/--use-raw are ignored and no gene scoring runs.")
    p.add_argument("--cluster-order", nargs="*", default=None)
    p.add_argument("--cluster-colors", default=None,
                   help="Optional JSON mapping cluster id -> hex colour (e.g. "
                        "data/cluster_colors_18.json). Colours the violin bodies and "
                        "x tick labels per cluster; unmapped clusters fall back to grey.")
    p.add_argument("--use-raw", action="store_true", default=True,
                   help="Score from .raw (lognorm). Default True — this h5ad has "
                        "scaled .X which produces all-NaN scores.")
    p.add_argument("--plot-format", nargs="+", default=["pdf", "png"])
    p.add_argument("--figwidth-mm", type=float, default=None,
                   help="author panel at this FINAL width in mm (assemble_figure --qa onpage_mm target) so it scales ~1:1; the saved PDF is wider by any outside legend — tune with --qa")
    p.add_argument("--figheight-mm", type=float, default=None, help="panel height in mm")
    return p.parse_args()


def load_gene_list(path):
    if path is None:
        return list(DEFAULT_PROLIF)
    with open(path) as fh:
        return [ln.strip() for ln in fh
                if ln.strip() and not ln.strip().startswith("#")]


def resolve_cluster_order(values):
    """Turn the --cluster-order argument into a flat list of string labels.

    Accepts either literal labels (e.g. `--cluster-order 11 15 2 ...`) or a single
    path to a JSON file. The JSON may be a bare list, or a dict with an "order" key
    (the canonical cluster_order.json layout). Labels are coerced to str.
    """
    if len(values) == 1 and values[0].lower().endswith(".json"):
        path = values[0]
        if not os.path.exists(path):
            raise SystemExit(f"--cluster-order JSON '{path}' not found")
        with open(path) as fh:
            obj = json.load(fh)
        order = obj["order"] if isinstance(obj, dict) else obj
        return [str(c) for c in order]
    return [str(c) for c in values]


def load_colors(path):
    """Flat {cluster_id: hex} map; keys coerced to str. Returns {} if unset/missing."""
    if not path or not os.path.exists(path):
        if path:
            print(f"  [warn] --cluster-colors '{path}' not found; using default colour")
        return {}
    with open(path) as fh:
        return {str(k): v for k, v in json.load(fh).items()}


def main():
    args = parse_args()
    plt = setup_style(font_size=8)
    os.makedirs(args.outdir, exist_ok=True)

    def figsize(default_w_in, default_h_in, w_mult=1.0, h_mult=1.0):
        fw, fh = args.figwidth_mm, args.figheight_mm
        if fw and fh:
            return ((fw / 25.4) * w_mult, (fh / 25.4) * h_mult)
        if fw:   # width only -> preserve authored aspect ratio (avoids skew + empty margins)
            w = (fw / 25.4) * w_mult
            return (w, default_h_in * w / default_w_in)
        if fh:   # height only -> preserve authored aspect ratio
            h = (fh / 25.4) * h_mult
            return (default_w_in * h / default_h_in, h)
        return (default_w_in, default_h_in)

    print(f"reading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    if args.cluster_col not in adata.obs.columns:
        raise SystemExit(
            f"--cluster-col '{args.cluster_col}' not in obs. Available: "
            + ", ".join(sorted(adata.obs.columns))
        )

    use_precomputed = bool(args.score_key)
    if use_precomputed:
        if args.score_key not in adata.obs.columns:
            raise SystemExit(
                f"--score-key '{args.score_key}' not in obs. Available: "
                + ", ".join(sorted(adata.obs.columns))
            )
        score_col = args.score_key
        print(f"using precomputed per-cell score obs['{score_col}'] "
              "(no gene scoring)")
    else:
        score_col = args.score_name
        genes_in = load_gene_list(args.gene_list)
        var_names = (adata.raw.var_names if (args.use_raw and adata.raw is not None)
                     else adata.var_names)
        present = [g for g in genes_in if g in var_names]
        missing = [g for g in genes_in if g not in var_names]
        print(f"proliferation genes: {len(present)}/{len(genes_in)} present, "
              f"{len(missing)} missing")
        if missing:
            print("  missing:", ", ".join(missing[:20])
                  + (" ..." if len(missing) > 20 else ""))

        used_df = pd.DataFrame({
            "gene": genes_in,
            "in_adata": [g in var_names for g in genes_in],
        })
        used_out = os.path.join(args.outdir, "fig3f_proliferation_genes_used.tsv")
        used_df.to_csv(used_out, sep="\t", index=False)
        print(f"  wrote {used_out}")

    data_labels = sorted(adata.obs[args.cluster_col].astype(str).unique(), key=numkey)
    if args.cluster_order is not None:
        requested = resolve_cluster_order(args.cluster_order)
        order = [c for c in requested if c in set(data_labels)]
        unknown = [c for c in requested if c not in set(data_labels)]
        dropped = [c for c in data_labels if c not in set(requested)]
        if unknown:
            print(f"  [warn] --cluster-order labels not in obs['{args.cluster_col}'] "
                  f"(ignored): {unknown}")
        if dropped:
            print(f"  [warn] obs labels missing from --cluster-order "
                  f"(not plotted): {dropped}")
        if not order:
            raise SystemExit(
                f"ERROR: none of the --cluster-order labels match obs"
                f"['{args.cluster_col}'].\n"
                f"  requested: {requested}\n"
                f"  actual labels in data: {data_labels}\n"
                "Pass labels from the same namespace (e.g. raw integer cluster ids), "
                "or drop --cluster-order to auto-order."
            )
    else:
        order = data_labels
    # categories MUST cover every label present, else pd.Categorical silently maps the
    # rest to NaN -> "nan" and they match no cluster. Append any extras after `order`.
    categories = order + [c for c in data_labels if c not in set(order)]
    adata.obs[args.cluster_col] = pd.Categorical(
        adata.obs[args.cluster_col].astype(str),
        categories=categories, ordered=True,
    )

    if not use_precomputed:
        print(f"sc.tl.score_genes -> obs['{score_col}'] (use_raw={args.use_raw})")
        sc.tl.score_genes(
            adata, gene_list=present, score_name=score_col,
            use_raw=args.use_raw, random_state=0,
        )
    nan_frac = adata.obs[score_col].isna().mean()
    if nan_frac > 0.5:
        raise SystemExit(
            f"ERROR: {nan_frac:.1%} of cells have NaN scores in obs['{score_col}']. "
            + ("The precomputed score column is mostly empty — check --score-key."
               if use_precomputed else
               "This usually means score_genes ran on scaled .X instead of lognorm .raw. "
               "Pass --use-raw (now the default) or check that .raw is populated.")
        )

    # Use numpy indexing instead of groupby — score_genes rewrites adata.obs internally
    # in some AnnData/scanpy versions, which can silently reset the Categorical dtype and
    # cause groupby(observed=True) to return no groups.
    _cl_arr = adata.obs[args.cluster_col].astype(str).values
    _sc_arr = adata.obs[score_col].values.astype(float)

    rows = []
    for _cl in order:
        _mask = _cl_arr == _cl
        _arr = _sc_arr[_mask]
        _valid = _arr[~np.isnan(_arr)]
        rows.append({
            "cluster": _cl,
            "count": int(_valid.size),
            "mean":   float(np.mean(_valid))           if _valid.size else float("nan"),
            "median": float(np.median(_valid))         if _valid.size else float("nan"),
            "q25":    float(np.percentile(_valid, 25)) if _valid.size else float("nan"),
            "q75":    float(np.percentile(_valid, 75)) if _valid.size else float("nan"),
        })
    summary = pd.DataFrame(rows)
    summary_out = os.path.join(args.outdir, "fig3f_proliferation_score_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"  wrote {summary_out}")

    colors = load_colors(args.cluster_colors)
    default_color = "#bdbdbd"

    by_cl = {r["cluster"]: _sc_arr[_cl_arr == r["cluster"]]
             for _, r in summary.iterrows() if r["count"] > 0}
    by_cl = {cl: arr[~np.isnan(arr)] for cl, arr in by_cl.items()}
    plot_order = [cl for cl in order if by_cl.get(cl, np.array([])).size > 0]
    vals = [by_cl[cl] for cl in plot_order]
    positions = list(range(len(plot_order)))

    # diagnostic: show exact sizes reaching violinplot
    print("  [debug] cluster score counts: " +
          " ".join(f"{cl}:{arr.size}" for cl, arr in zip(plot_order, vals)))

    # hard guard: violinplot crashes on arrays with <2 elements in some mpl versions
    _ok = [(cl, arr) for cl, arr in zip(plot_order, vals) if arr.size >= 2]
    if not _ok:
        print("  WARNING: no clusters have >=2 non-NaN scores; skipping violin plot")
        return
    _dropped = [cl for cl, arr in zip(plot_order, vals) if arr.size < 2]
    if _dropped:
        print(f"  [warn] skipped cluster(s) with <2 scores: {_dropped}")
    plot_order, vals = map(list, zip(*_ok))
    positions = list(range(len(plot_order)))

    fig, ax = plt.subplots(figsize=figsize(max(5.0, len(plot_order) * 0.4), 3.2))
    parts = ax.violinplot(vals, positions=positions, widths=0.8,
                          showmeans=False, showextrema=False, showmedians=True)
    for body, cl in zip(parts["bodies"], plot_order):
        body.set_facecolor(colors.get(cl, default_color))
        body.set_edgecolor("black")
        body.set_linewidth(0.3)
        body.set_alpha(0.9)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    ax.axhline(0, color="#888", lw=0.5, ls="--")
    ax.set_xticks(positions)
    # cluster IDs vertical: short labels stay readable and never collide at narrow widths
    ax.set_xticklabels(plot_order, rotation=90, fontsize=7)
    if colors:
        for t, cl in zip(ax.get_xticklabels(), plot_order):
            t.set_color(colors.get(cl, "black"))
    ax.set_ylabel("proliferation score")
    ax.set_xlabel("cluster")# (adenocarcinoma -> neuroendocrine)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()

    savefig(fig, args.outdir, "fig3f_proliferation_violin",
            formats=args.plot_format)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
