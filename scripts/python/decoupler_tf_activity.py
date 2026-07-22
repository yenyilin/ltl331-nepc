#!/usr/bin/env python3
"""
decoupler_tf_activity.py — orthogonal prior-knowledge TF-activity layer
(decoupleR + CollecTRI/DoRothEA) for the LTL331 trajectory. Cross-validates
SCENIC findings and checks whether canonical NEPC TFs (ASCL1, FOXA2, NKX2-1,
POU3F2/BRN2, NEUROD1, ONECUT2, INSM1) score as active in cluster 17 / NE
clusters even when not in the 200-regulon SCENIC set.

Method
  1. Build per-group pseudobulk = mean lognorm expression. Default group =
     seurat_clusters; pass --pseudobulk-by 'seurat_clusters,sample' for the
     per-cluster x per-timepoint dynamics.
  2. Load prior network (CollecTRI by default via OmniPath, or DoRothEA; falls
     back to a local --prior-tsv if the analysis node is air-gapped).
  3. Pre-filter expression to genes in the prior (huge memory win); drop TFs
     whose regulon doesn't overlap the expression matrix by >= --min-target-overlap.
  4. Run ULM (default; MLM and WSUM also exposed).
  5. Emit TF x cluster activity TSV, heatmap of top-variance TFs, per-cluster
     top-K TFs, and a focused canonical-NEPC-TF readout for the response.

Inputs
  --adata        velocity h5ad (seurat_clusters + cellid + sample)
  --full-h5ad    full-gene h5ad with lognorm in .raw (or .X)

Prior options
  --prior collectri     (default; fetched from OmniPath; needs internet once)
  --prior dorothea      (older; confidence levels A,B,C kept)
  --prior-tsv X.tsv     (load from local TSV — bypasses network call)

Outputs (in --outdir)
  tf_activity_per_cluster.tsv      TF × cluster activity (ULM/MLM estimate)
  tf_activity_pvals.tsv            TF × cluster p-values (when method returns them)
  tf_activity_heatmap.{pdf,png}    Top-N most-variable TFs
  top_tfs_per_cluster.tsv          Top-K TFs by |activity| per cluster
  canonical_nepc_tfs.tsv           Focused readout: ASCL1/FOXA2/NKX2-1/POU3F2/...
  prior_used.tsv                   The exact filtered prior used (provenance)

Example
  python3 scripts/decoupler_tf_activity.py \\
      --adata ltl331_velocity_dynamical.h5ad \\
      --full-h5ad full_h5ad.h5ad \\
      --pseudobulk-by seurat_clusters \\
      --prior collectri --method ulm \\
      --top-n 30 \\
      --outdir decoupler_tf
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from _utils import setup_style, savefig


# canonical TFs we want a focused readout for, regardless of overall ranking
CANONICAL_NEPC_TFS = [
    # NEPC drivers
    'ASCL1', 'NEUROD1', 'POU3F2', 'POU2F3', 'INSM1', 'FOXA2', 'NKX2-1',
    'ONECUT2', 'REST', 'SOX2', 'MSX1', 'YAP1',
    # PRAD / castration-resistance reference points
    'AR', 'FOXA1', 'NKX3-1', 'TP63', 'TP53', 'MYC',
]


def get_lognorm(full):
    """Return the lognorm-bearing AnnData (full.raw -> adata if present)."""
    if full.raw is not None:
        try:
            return full.raw.to_adata()
        except Exception:
            pass
    return full


def to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def pseudobulk_mean(adata_full, cells_in_order, groups):
    """Mean log-expression per group. Returns (groups × genes) DataFrame."""
    sub = adata_full[cells_in_order, :]
    X = to_dense(sub.X)
    df = pd.DataFrame(X, columns=sub.var_names.tolist())
    df['_grp'] = pd.Series(groups, dtype='category').values
    return df.groupby('_grp').mean()


def load_prior(args):
    """Load TF -> target prior network as a DataFrame with [source,target,weight].

    Uses decoupler-py 2.x API (`dc.op.collectri` / `dc.op.dorothea`), with a
    1.x fallback (`dc.get_collectri` / `dc.get_dorothea`) for older installs.
    Verified against decoupler 2.1.6 source.
    """
    if args.prior_tsv:
        sep = '\t' if str(args.prior_tsv).endswith(('.tsv', '.txt')) else ','
        net = pd.read_csv(args.prior_tsv, sep=sep)
    else:
        try:
            import decoupler as dc
        except ImportError:
            raise SystemExit(
                'decoupler not importable. Install: pip install decoupler '
                '(PyPI distribution is "decoupler", NOT "decoupler-py")')
        is_v2 = hasattr(dc, 'op')
        ver = getattr(dc, '__version__', '?')
        print(f'decoupler API: {"2.x" if is_v2 else "1.x"} (v{ver})')
        if args.prior == 'collectri':
            print(f'fetching CollecTRI (license={args.license})')
            if is_v2:
                # 2.x signature: collectri(organism, remove_complexes, license, verbose)
                net = dc.op.collectri(organism='human',
                                      remove_complexes=False,
                                      license=args.license)
            else:
                # 1.x signature: get_collectri(organism, split_complexes)
                net = dc.get_collectri(organism='human', split_complexes=False)
        elif args.prior == 'dorothea':
            print(f'fetching DoRothEA (levels A,B,C; license={args.license})')
            if is_v2:
                net = dc.op.dorothea(organism='human',
                                     levels=['A', 'B', 'C'],
                                     license=args.license)
            else:
                net = dc.get_dorothea(levels=['A', 'B', 'C'], organism='human')
        else:
            raise SystemExit(f'unknown --prior {args.prior}')

    # normalize column names (older fetchers used 'mor' / 'tf')
    rename = {}
    if 'mor' in net.columns and 'weight' not in net.columns:
        rename['mor'] = 'weight'
    if 'sign' in net.columns and 'weight' not in net.columns:
        rename['sign'] = 'weight'
    if 'tf' in net.columns and 'source' not in net.columns:
        rename['tf'] = 'source'
    net = net.rename(columns=rename)
    if 'weight' not in net.columns:
        print('[note] no weight column in prior; defaulting to +1 (treating all '
              'interactions as activating)')
        net['weight'] = 1
    needed = {'source', 'target', 'weight'}
    if not needed.issubset(net.columns):
        raise SystemExit(f'prior missing required columns: have {list(net.columns)}, '
                         f'need {needed}')
    print(f'prior network loaded: {len(net):,} interactions; '
          f'{net["source"].nunique():,} TFs; {net["target"].nunique():,} target genes')
    return net


def run_method(method, mat, net, tmin=1):
    """Run decoupler activity method on a pseudobulk DataFrame.

    Uses decoupler-py 2.x API (`dc.mt.<method>`) by default; 1.x fallback
    (`dc.run_<method>`) for older installs. Verified against decoupler 2.1.6.

    In 2.x:
      - Signature: method(data, net, tmin=5, raw=False, empty=True,
                          bsize=250_000, verbose=False, **kwargs)
        First positional is `data`; net's source/target/weight columns are
        auto-detected by decoupler — no need to pass them.
      - `wsum` was renamed to `waggr(fun="wsum")`; `wmean` is `waggr(fun="wmean")`.
      - For DataFrame input, returns `(es, pv)` tuple.

    `tmin=1` here because we already pre-filter the prior in main(); pass a
    larger value to let decoupler also enforce a regulon-size floor.

    Returns (estimate_df, pvals_df_or_None).
    """
    import decoupler as dc
    is_v2 = hasattr(dc, 'mt')
    if is_v2:
        if method == 'wsum':
            fn = dc.mt.waggr
            extra = {'fun': 'wsum'}
            print('[note] in decoupler 2.x, wsum is dc.mt.waggr(fun="wsum")')
        elif method == 'wmean':
            fn = dc.mt.waggr
            extra = {'fun': 'wmean'}
        else:
            fn = getattr(dc.mt, method, None)
            extra = {}
        if fn is None:
            raise SystemExit(
                f'method "{method}" not in decoupler.mt. Available: '
                f'ulm, mlm, zscore, viper, aucell, gsea, gsva, mdt, udt, ora, '
                f'waggr (use wsum/wmean via --method).'
            )
        out = fn(data=mat, net=net, tmin=tmin, verbose=False, **extra)
    else:
        fn = getattr(dc, f'run_{method}', None)
        if fn is None:
            raise SystemExit(
                f'method "{method}" not available in decoupler 1.x. '
                f'Try ulm, mlm, or wsum.'
            )
        out = fn(mat=mat, net=net, source='source', target='target',
                 weight='weight', verbose=False)

    # DataFrame input -> tuple (es, pv) in both 2.x and 1.x
    if isinstance(out, tuple):
        if len(out) >= 2:
            return out[0], out[1]
        if len(out) == 1:
            return out[0], None
        return None, None
    if isinstance(out, pd.DataFrame):
        return out, None
    # AnnData write-back (only if caller passed an AnnData)
    if out is None and hasattr(mat, 'obsm'):
        name = method if method not in ('wsum', 'wmean') else 'waggr'
        sk = f'score_{name}'
        if sk in mat.obsm:
            est = pd.DataFrame(mat.obsm[sk], index=mat.obs_names)
            pv = (pd.DataFrame(mat.obsm[f'padj_{name}'], index=mat.obs_names)
                  if f'padj_{name}' in mat.obsm else None)
            return est, pv
    raise SystemExit(f'{method}: unexpected output type {type(out).__name__}')


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--adata', required=True)
    p.add_argument('--full-h5ad', required=True)
    p.add_argument('--cellid-key', default='cellid')
    p.add_argument('--min-match-frac', type=float, default=0.5,
                   help='abort with a barcode-convention diagnostic if fewer than this '
                        'fraction of --cellid-key values match --full-h5ad obs_names '
                        '(guards the historical 0/51,726 silent-mismatch crash)')
    p.add_argument('--pseudobulk-by', default='seurat_clusters',
                   help='comma-separated obs cols to group by; e.g. seurat_clusters '
                        'or "seurat_clusters,sample" for cluster × timepoint dynamics')
    p.add_argument('--prior', choices=['collectri', 'dorothea'], default='collectri')
    p.add_argument('--prior-tsv', default=None,
                   help='load prior from local TSV instead of OmniPath '
                        '(columns: source,target,weight [,confidence])')
    p.add_argument('--license', default='academic',
                   choices=['academic', 'nonprofit', 'commercial'],
                   help='OmniPath license tier for decoupler 2.x prior fetchers '
                        '(ignored on 1.x); default: academic')
    p.add_argument('--method', choices=['ulm', 'mlm', 'wsum'], default='ulm')
    p.add_argument('--min-target-overlap', type=int, default=5,
                   help='drop TFs whose regulon overlaps the expression matrix '
                        'by < this many genes')
    p.add_argument('--top-n', type=int, default=30,
                   help='top N TFs (by activity variance across groups) for the heatmap')
    p.add_argument('--top-k-per-cluster', type=int, default=10)
    p.add_argument('--canonical-tfs', nargs='+', default=CANONICAL_NEPC_TFS,
                   help='TFs to include in the focused cross-validation table')
    p.add_argument('--zscore-pseudobulk', action='store_true',
                   help='z-score each gene across pseudobulk groups before running '
                        'decoupler methods. Removes the broad-regulon bias that '
                        'makes MYC/AR/TP53 universally positive and FOXA2 etc. '
                        'spuriously negative on raw lognorm pseudobulks. STRONGLY '
                        'recommended for cross-cluster TF-activity comparisons.')
    p.add_argument('--all-methods', action='store_true',
                   help='also run the other two methods (ULM + MLM + WSUM) and '
                        'produce a method-concordance table for --canonical-tfs')
    p.add_argument('--robust-rank-cutoff', type=int, default=50,
                   help='under --all-methods, a canonical TF × cluster pair is '
                        'flagged "robust" if its |activity| rank within the cluster '
                        'is <= this value under ALL methods (default 50)')
    p.add_argument('--outdir', default='decoupler_tf')
    p.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # 1. load data
    import scanpy as sc
    adata = sc.read_h5ad(args.adata)
    full = sc.read_h5ad(args.full_h5ad)
    full_ln = get_lognorm(full)
    xmin = to_dense(full_ln[:5].X).min()
    if xmin < -0.01:
        print(f'[warn] full source X looks scaled (min={xmin:.2f}); ULM is '
              f'calibrated for log-norm input. Rankings still meaningful but '
              f'absolute activity values are on a scaled scale.')

    # 2. align cellid → full_ln
    key = adata.obs[args.cellid_key].astype(str)
    in_full = key.isin(full_ln.obs_names).values
    frac = in_full.mean() if len(key) else 0.0
    print(f'cellid-matched cells: {in_full.sum():,}/{len(key):,} ({frac:.1%})')
    # fail-fast on a barcode-CONVENTION mismatch (the historical 0/51,726 bug) rather
    # than crashing 200 lines later in the groupby with a cryptic Categorical error.
    if frac < args.min_match_frac:
        import re as _re
        def _bare(s):
            m = _re.search(r'[ACGTacgt]{16}', str(s)); return m.group(0).upper() if m else '?'
        ex_key = key.head(3).tolist()
        ex_full = list(full_ln.obs_names[:3])
        # quick auto-probe: would the bare 16-nt barcode match better?
        bare = key.map(_bare)
        full_bare = {_bare(x) for x in full_ln.obs_names[:50000]}
        hint = (f'  bare-16nt barcodes would match ~{bare.isin(full_bare).mean():.0%} '
                f'— try canonicalize_cellid.py to reconcile conventions'
                if bare.isin(full_bare).mean() > frac else
                '  even bare barcodes do not match — check you passed the right objects/genome')
        raise SystemExit(
            f'\n[abort] only {frac:.1%} of --cellid-key ({args.cellid_key!r}) values match '
            f'--full-h5ad obs_names (threshold --min-match-frac={args.min_match_frac}).\n'
            f'  This is a barcode-CONVENTION mismatch, not a biology problem.\n'
            f'  adata.obs[{args.cellid_key!r}][:3] = {ex_key}\n'
            f'  full_ln.obs_names[:3]            = {ex_full}\n'
            f'{hint}\n'
            f'  Fix: run  scripts/canonicalize_cellid.py diagnose --adata <velocity> '
            f'--target <full> , then `apply`, OR pass a different --cellid-key.\n')
    cells = key[in_full].tolist()

    # 3. group labels (single column or composite)
    by_cols = [c.strip() for c in args.pseudobulk_by.split(',')]
    for c in by_cols:
        if c not in adata.obs.columns:
            raise SystemExit(f'--pseudobulk-by column not in adata.obs: {c}')
    groups = (adata.obs.loc[in_full, by_cols].astype(str)
              .agg('|'.join, axis=1).values)
    n_groups = len(np.unique(groups))
    print(f'pseudobulks: grouping by {by_cols} -> {n_groups} group(s)')

    # 4. load prior, then SUBSET expression matrix to genes in the prior
    #    (big memory + speed win; ULM only sees those genes anyway)
    net = load_prior(args)
    target_set = set(net['target'].unique())
    keep_genes = [g for g in full_ln.var_names if g in target_set]
    print(f'expression genes also in prior targets: {len(keep_genes):,} '
          f'(out of {full_ln.n_vars:,} full)')
    if len(keep_genes) < 100:
        print('[warn] very few genes shared between expression and prior — '
              'check organism (CollecTRI is human) and gene symbol convention')
    full_sub = full_ln[:, keep_genes]

    # 5. pseudobulk on the subsetted matrix
    print(f'building pseudobulk: {len(cells):,} cells × {len(keep_genes):,} '
          f'genes (densified for groupby)')
    pdata = pseudobulk_mean(full_sub, cells, groups)
    pdata.index.name = '|'.join(by_cols)
    print(f'pseudobulk matrix: {pdata.shape[0]} groups × {pdata.shape[1]} genes')

    # 5b. optional: z-score each gene across clusters → makes ULM scores reflect
    #     CLUSTER-LEVEL DEVIATION rather than absolute pseudobulk magnitude. This
    #     removes the broad-regulon bias (MYC/AR/TP53 universally large, FOXA2 etc.
    #     spuriously negative) and is the standard fix recommended in the decoupler
    #     tutorials for cross-cluster TF-activity comparisons.
    if args.zscore_pseudobulk:
        print('z-scoring pseudobulk per gene across groups (cross-cluster '
              'deviation scale; recommended for ULM/MLM interpretation)')
        mu = pdata.mean(axis=0)
        sd = pdata.std(axis=0).replace(0, 1)
        pdata = pdata.sub(mu, axis=1).div(sd, axis=1)
        print(f'  post-zscore range: min={pdata.values.min():+.2f}  '
              f'max={pdata.values.max():+.2f}  '
              f'mean={pdata.values.mean():+.4f} (≈0 expected)')
    else:
        print('[note] pseudobulk NOT z-scored. ULM/MLM activity scores will be '
              'biased by absolute baseline expression and regulon size '
              '(broad-regulon TFs like MYC/AR/TP53 will be uniformly positive; '
              'small or cross-tissue regulons like FOXA2 may be uniformly '
              'negative). Interpret BY WITHIN-TF RANKING across clusters. To '
              'get cleaner per-cluster activity scores, re-run with '
              '--zscore-pseudobulk.')

    # 6. filter prior: keep only TFs with >= min_target_overlap targets in expr
    avail = set(pdata.columns)
    net = net[net['target'].isin(avail)].copy()
    tf_target_n = net.groupby('source').size()
    keep_tfs = tf_target_n[tf_target_n >= args.min_target_overlap].index
    net = net[net['source'].isin(keep_tfs)].copy()
    print(f'prior after TF/target filtering: '
          f'{net["source"].nunique():,} TFs, {len(net):,} interactions')
    net.to_csv(outdir / 'prior_used.tsv', sep='\t', index=False)

    # 7. run decoupler
    print(f'\nrunning {args.method.upper()} ...')
    est, pv = run_method(args.method, pdata, net)
    print(f'activity matrix: {est.shape[0]} groups × {est.shape[1]} TFs')

    est.to_csv(outdir / 'tf_activity_per_cluster.tsv', sep='\t')
    print(f'wrote {outdir}/tf_activity_per_cluster.tsv')
    if pv is not None:
        pv.to_csv(outdir / 'tf_activity_pvals.tsv', sep='\t')
        print(f'wrote {outdir}/tf_activity_pvals.tsv')

    # 8. top-N variance TFs for heatmap
    var = est.var(axis=0).sort_values(ascending=False)
    top_tfs = var.head(args.top_n).index.tolist()
    print(f'\ntop-{args.top_n} most variable TFs (first 10): {top_tfs[:10]}')

    # 9. heatmap (TFs × groups)
    plt = setup_style(font_size=7)
    H = est[top_tfs].T
    fig, ax = plt.subplots(figsize=(0.5 * H.shape[1] + 2, 0.20 * H.shape[0] + 2))
    vmax = float(np.nanpercentile(np.abs(H.values), 95))
    im = ax.imshow(H.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    ax.set_xticks(range(H.shape[1]))
    ax.set_xticklabels(H.columns, rotation=45, ha='right')
    ax.set_yticks(range(H.shape[0]))
    ax.set_yticklabels(H.index)
    plt.colorbar(im, ax=ax, fraction=0.03,
                 label=f'{args.method.upper()} activity ({args.prior})')
    ax.set_title(f'Top-{args.top_n} variable TFs (decoupleR / {args.prior} / '
                 f'{args.method.upper()})')
    savefig(fig, outdir, 'tf_activity_heatmap', args.plot_format)
    plt.close(fig)

    # 10. top-K TFs per group (by |activity|)
    rows = []
    for grp in est.index:
        top = est.loc[grp].reindex(est.loc[grp].abs().sort_values(ascending=False)
                                    .head(args.top_k_per_cluster).index)
        for rank, (tf, act) in enumerate(top.items(), 1):
            rows.append({
                'group': grp, 'rank': rank, 'tf': tf,
                'activity': float(act),
                'pval': (float(pv.loc[grp, tf])
                         if pv is not None and tf in pv.columns else np.nan),
            })
    pd.DataFrame(rows).to_csv(outdir / 'top_tfs_per_cluster.tsv', sep='\t', index=False)
    print(f'wrote {outdir}/top_tfs_per_cluster.tsv')

    # 11. canonical NEPC TFs focused readout
    cano_present = [t for t in args.canonical_tfs if t in est.columns]
    cano_missing = [t for t in args.canonical_tfs if t not in est.columns]
    if cano_missing:
        print(f'\n[note] canonical TFs not scored (no regulon in filtered prior, '
              f'or insufficient target overlap): {cano_missing}')
    if cano_present:
        cano_df = est[cano_present].T
        cano_df.index.name = 'TF'
        cano_df.to_csv(outdir / 'canonical_nepc_tfs.tsv', sep='\t')
        print(f'\ncanonical NEPC TF activity (rows = TF, cols = group):')
        print(cano_df.round(2).to_string())

    # 12. optional method-concordance check (--all-methods)
    if args.all_methods:
        print('\n=== --all-methods: running ULM + MLM + WSUM for concordance ===')
        from scipy.stats import spearmanr
        all_est = {args.method: est}
        for m in ['ulm', 'mlm', 'wsum']:
            if m == args.method:
                continue
            try:
                print(f'\nrunning {m.upper()} ...')
                e_m, pv_m = run_method(m, pdata, net)
                e_m.to_csv(outdir / f'tf_activity_{m}.tsv', sep='\t')
                if pv_m is not None:
                    pv_m.to_csv(outdir / f'tf_activity_{m}_pvals.tsv', sep='\t')
                print(f'  -> {e_m.shape[0]} groups × {e_m.shape[1]} TFs; '
                      f'wrote tf_activity_{m}.tsv')
                all_est[m] = e_m
            except Exception as exc:
                print(f'[error] {m.upper()} failed: {exc}; skipping it for concordance')
        # also rename primary output for consistency
        est.to_csv(outdir / f'tf_activity_{args.method}.tsv', sep='\t')

        # method-pair Spearman correlation across all (group × TF) entries
        methods_ok = list(all_est.keys())
        pair_rows = []
        for i, m1 in enumerate(methods_ok):
            for m2 in methods_ok[i + 1:]:
                e1, e2 = all_est[m1], all_est[m2]
                common_tfs = e1.columns.intersection(e2.columns)
                common_grp = e1.index.intersection(e2.index)
                v1 = e1.loc[common_grp, common_tfs].values.flatten()
                v2 = e2.loc[common_grp, common_tfs].values.flatten()
                ok = np.isfinite(v1) & np.isfinite(v2)
                rho, p_ = spearmanr(v1[ok], v2[ok])
                pair_rows.append({
                    'method_a': m1, 'method_b': m2,
                    'n_TFs_common': int(len(common_tfs)),
                    'n_entries': int(ok.sum()),
                    'spearman_rho': float(rho), 'spearman_p': float(p_),
                })
                print(f'concordance {m1.upper()} vs {m2.upper()}: '
                      f'rho={rho:+.3f}  (n={int(ok.sum())} TF×group pairs)')
        pd.DataFrame(pair_rows).to_csv(outdir / 'method_pair_correlation.tsv',
                                       sep='\t', index=False)

        # per-canonical-TF × cluster concordance: rank |activity| per cluster,
        # flag pair as "robust" if rank <= robust_rank_cutoff under ALL methods
        K = args.robust_rank_cutoff
        conc_rows = []
        for tf in args.canonical_tfs:
            for grp in pdata.index:
                row = {'TF': tf, 'cluster': grp}
                ranks = []
                for m in methods_ok:
                    e = all_est[m]
                    if tf not in e.columns or grp not in e.index:
                        row[f'{m}_activity'] = np.nan
                        row[f'{m}_rank_abs'] = np.nan
                        ranks.append(np.nan)
                        continue
                    a = e.loc[grp, tf]
                    row[f'{m}_activity'] = float(a)
                    abs_act = e.loc[grp].abs()
                    rank = int((abs_act > abs(a)).sum() + 1)
                    row[f'{m}_rank_abs'] = rank
                    ranks.append(rank)
                valid = [r for r in ranks if not (isinstance(r, float) and np.isnan(r))]
                row['robust_top_k'] = (
                    len(valid) == len(methods_ok) and all(r <= K for r in valid)
                )
                conc_rows.append(row)
        conc_df = pd.DataFrame(conc_rows)
        conc_df.to_csv(outdir / 'method_concordance_canonical.tsv',
                       sep='\t', index=False)

        n_robust = int(conc_df['robust_top_k'].sum())
        n_pairs = len(conc_df)
        print(f'\ncanonical TF × cluster pairs ROBUST across all methods '
              f'(top-{K} by |activity|): {n_robust} / {n_pairs}')
        # by-TF summary
        by_tf = (conc_df.groupby('TF')['robust_top_k'].sum()
                 .sort_values(ascending=False))
        print('clusters with robust call, per canonical TF:')
        print(by_tf.to_string())

        print(f'\nwrote {outdir}/method_pair_correlation.tsv')
        print(f'wrote {outdir}/method_concordance_canonical.tsv')
        for m in methods_ok:
            print(f'wrote {outdir}/tf_activity_{m}.tsv')

    print('\nReading guide:')
    print('  tf_activity_per_cluster.tsv  full TF × cluster activity (signed; primary method).')
    print(f'  tf_activity_heatmap.*        top-{args.top_n} most variable TFs (primary method).')
    print('  top_tfs_per_cluster.tsv      top-K TFs per cluster for the response.')
    print('  canonical_nepc_tfs.tsv       ASCL1/FOXA2/NKX2-1/POU3F2/NEUROD1/ONECUT2')
    print('                                in cluster 17 / NE clusters — directly')
    print('                                answers the TF-related question without re-running')
    print('                                SCENIC.')
    if args.all_methods:
        print('  tf_activity_<method>.tsv     activity matrices per method (ulm/mlm/wsum).')
        print('  method_pair_correlation.tsv  pairwise Spearman concordance across methods.')
        print('  method_concordance_canonical.tsv')
        print(f'                                per-canonical-TF × cluster: activity + rank under')
        print(f'                                each method + "robust_top_k" flag (rank ≤ {args.robust_rank_cutoff} in all).')
    print('  Cross-validation hint: TFs that score in the top-N here AND appear as')
    print('  top SCENIC regulons (or vice versa) are the high-confidence drivers.')


if __name__ == '__main__':
    main()
