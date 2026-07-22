#!/usr/bin/env python3
"""
clinical_signature_validation.py — LTL331→clinical signature validation at
the PATIENT level (not per-dataset clusters, which were unstable).

What it does
------------
1. Loads LTL331 per-cluster signatures (from make_signatures.py: cluster -> gene list).
2. Scores every signature on every cell with scanpy `sc.tl.score_genes` on .raw lognorm
   (the Python equivalent of Seurat AddModuleScore/UCell so no R).
3. Aggregates to PATIENT level (mean score per clinical patient), keyed by orig.ident,
   and joins the verified disease phenotype (patient_phenotype / ned_status) + a disease
   axis ordinal (PIN < HSPC < PRAD < CRPC < NEPC by default; overridable).
4. Tests each LTL331 signature for a monotonic trend along the disease axis
   (Spearman + Jonckheere–Terpstra with a permutation p), and NED vs non-NED
   (Mann–Whitney). n~10 patients -> effect size first, p as support.
5. Emits the disease-axis-ordered patient x signature heatmap (z-scored per signature)
   + all stats as TSV. This replaces the old per-cluster violins.

Why patient-level: per-dataset clinical clusters (Gao_0/Gao_20/…) are not cross-cohort
comparable and shift whenever clustering is redone. Patient + verified pathology is stable.

Deps: scanpy, pandas, numpy, scipy, matplotlib (all Python). No R, no UCell.

Usage
  python clinical_signature_validation.py \\
      --h5ad harmony.annotated.h5ad --signatures ltl331_cluster_markers.tsv \\
      --out data/clinical_validation --group-col dataset --ltl331-label LTL331 \\
      --patient-col orig.ident --phenotype-col patient_phenotype \\
      --cluster-order data/cluster_order.json --plot-format pdf png
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import setup_style, savefig, numkey  # noqa: E402

# default patient_phenotype -> (disease_state, disease_axis ordinal). Override with --metadata.
DEFAULT_AXIS = {
    'Gao_PIN':                ('PIN',  0),
    'Li_mHSPC_adeno':         ('HSPC', 1),
    'Li_nmHSPC_mixedNEPC':    ('HSPC', 1),
    'Gao_PRAD':               ('PRAD', 2),
    'Li_mCRPC_mixedNEPC':     ('CRPC', 3),
    'Gao_NEPC_smallcell':     ('NEPC', 4),
    'Li_mCRPC_pureSmallCell': ('NEPC', 4),
}


def load_signatures(path):
    """cluster -> [genes]; TSV lines '<cluster>\\t gene1,gene2,...'."""
    sigs = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cl, genes = line.rstrip('\n').split('\t', 1)
            sigs[str(cl)] = [g for g in genes.replace(',', ' ').split() if g]
    return sigs


def jonckheere(groups):
    """Jonckheere–Terpstra statistic for ordered groups (list of 1D arrays, in order)."""
    jt = 0.0
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a, b = groups[i], groups[j]
            # count pairs (x in a, y in b) with y > x, +0.5 ties
            jt += sum((np.sum(b > x) + 0.5 * np.sum(b == x)) for x in a)
    return jt


def jt_perm_p(values, axis, n_perm, rng):
    """Permutation p for a positive (increasing) JT trend."""
    order = np.argsort(axis, kind='stable')
    uniq = sorted(set(axis))
    def stat(v):
        groups = [v[np.array(axis) == u] for u in uniq]
        return jonckheere(groups)
    obs = stat(values)
    ge = 1
    for _ in range(n_perm):
        if stat(rng.permutation(values)) >= obs:
            ge += 1
    return obs, ge / (n_perm + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True, help='full integrated object (needs .raw lognorm)')
    ap.add_argument('--signatures', required=True, help='from make_signatures.py')
    ap.add_argument('--out', default='data/clinical_validation')
    ap.add_argument('--group-col', default='dataset')
    ap.add_argument('--ltl331-label', default='LTL331')
    ap.add_argument('--patient-col', default='orig.ident')
    ap.add_argument('--phenotype-col', default='patient_phenotype')
    ap.add_argument('--ned-col', default='ned_status')
    ap.add_argument('--cluster-order', default=None, help='order signature columns by this')
    ap.add_argument('--metadata', default=None,
                    help='optional TSV (patient_phenotype, disease_state, disease_axis) '
                         'overriding DEFAULT_AXIS')
    ap.add_argument('--use-raw', action='store_true', default=True)
    ap.add_argument('--no-raw', dest='use_raw', action='store_false')
    ap.add_argument('--n-perm', type=int, default=999)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--plot-format', nargs='+', default=['pdf', 'png'])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f'[load] {args.h5ad}')
    adata = sc.read_h5ad(args.h5ad)
    src = set(adata.raw.var_names if (args.use_raw and adata.raw is not None) else adata.var_names)

    # ---- score each signature (scanpy score_genes on .raw) -----------------
    sigs = load_signatures(args.signatures)
    order = sorted(sigs, key=numkey)
    if args.cluster_order and os.path.exists(args.cluster_order):
        co = [str(c) for c in json.load(open(args.cluster_order))['order']]
        pos = {c: i for i, c in enumerate(co)}
        order = sorted(sigs, key=lambda c: (pos.get(str(c), 1e9), numkey(c)))
    score_cols = []
    print('[score] sc.tl.score_genes on', '.raw' if args.use_raw else '.X')
    for cl in order:
        present = [g for g in sigs[cl] if g in src]
        if not present:
            print(f'   sig {cl}: 0 genes present — skipped'); continue
        name = f'sig_{cl}'
        sc.tl.score_genes(adata, present, score_name=name, use_raw=args.use_raw)
        score_cols.append(name)
        if len(present) < len(sigs[cl]):
            print(f'   sig {cl}: {len(present)}/{len(sigs[cl])} genes present')

    # ---- patient-level aggregation (clinical patients only) ----------------
    obs = adata.obs
    clin = obs[obs[args.group_col].astype(str) != args.ltl331_label].copy()
    pat = clin.groupby(args.patient_col, observed=True)[score_cols].mean()
    # attach phenotype (one per patient) + disease axis
    phen = clin.groupby(args.patient_col, observed=True)[args.phenotype_col].agg(
        lambda s: s.astype(str).mode().iat[0])
    ned = clin.groupby(args.patient_col, observed=True)[args.ned_col].agg(
        lambda s: s.astype(str).mode().iat[0]) if args.ned_col in clin else None
    axis_map = dict(DEFAULT_AXIS)
    if args.metadata and os.path.exists(args.metadata):
        m = pd.read_csv(args.metadata, sep=None, engine='python')
        axis_map.update({r[args.phenotype_col]: (r['disease_state'], int(r['disease_axis']))
                         for _, r in m.iterrows()})
    pat['phenotype'] = phen
    pat['disease_state'] = pat['phenotype'].map(lambda p: axis_map.get(p, ('?', 99))[0])
    pat['disease_axis'] = pat['phenotype'].map(lambda p: axis_map.get(p, ('?', 99))[1])
    if ned is not None:
        pat['ned_status'] = ned
    unmapped = pat.index[pat['disease_axis'] == 99].tolist()
    if unmapped:
        print(f'[warn] patients with unmapped phenotype (axis=99): '
              f'{pat.loc[unmapped, "phenotype"].to_dict()} — pass --metadata to place them')
    pat = pat.sort_values(['disease_axis', 'disease_state'])
    pat.to_csv(f'{args.out}/patient_signature_scores.tsv', sep='\t')
    print(f'[ok] {len(pat)} clinical patients x {len(score_cols)} signatures')

    # ---- trend tests per signature -----------------------------------------
    from scipy.stats import spearmanr, mannwhitneyu
    ax = pat['disease_axis'].to_numpy(float)
    rows = []
    for c in score_cols:
        v = pat[c].to_numpy(float)
        rho, p_s = spearmanr(v, ax)
        jt, p_jt = jt_perm_p(v, list(ax), args.n_perm, rng)
        row = {'signature': c.replace('sig_', ''), 'spearman_rho': round(float(rho), 3),
               'spearman_p': f'{p_s:.3g}', 'JT_stat': jt, 'JT_perm_p': f'{p_jt:.3g}'}
        if 'ned_status' in pat:
            nd = pat['ned_status'].astype(str)
            a = pat.loc[nd == 'NED', c]; b = pat.loc[nd == 'non-NED', c]
            if len(a) and len(b):
                u, p_u = mannwhitneyu(a, b, alternative='two-sided')
                row['NED_minus_nonNED'] = round(float(a.mean() - b.mean()), 3)
                row['NED_MWU_p'] = f'{p_u:.3g}'
        rows.append(row)
    trend = pd.DataFrame(rows)
    trend.to_csv(f'{args.out}/signature_trend_tests.tsv', sep='\t', index=False)
    print('\n[trend] per-signature monotonic trend along disease axis:')
    print(trend.to_string(index=False))

    # ---- heatmap: patients (rows, disease-ordered) x signatures (z-scored) -
    plt = setup_style(font_size=7)
    Z = pat[score_cols].copy()
    Z = (Z - Z.mean(0)) / Z.std(0).replace(0, np.nan)        # z per signature across patients
    fig_h = max(3.5, 0.32 * len(pat) + 2.0)
    fig_w = max(6.0, 0.42 * len(score_cols) + 3.0)
    fig, ax_h = plt.subplots(figsize=(fig_w, fig_h))
    im = ax_h.imshow(Z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2,
                     interpolation='nearest')
    ax_h.set_xticks(range(len(score_cols)))
    ax_h.set_xticklabels([c.replace('sig_', '') for c in score_cols], fontsize=9)
    ax_h.set_xlabel('LTL331 cluster signature (adenocarcinoma → neuroendocrine)', fontsize=11)
    ylab = [f'{p}  [{pat.loc[p, "disease_state"]}]' for p in pat.index]
    ax_h.set_yticks(range(len(pat))); ax_h.set_yticklabels(ylab, fontsize=9)
    ax_h.set_ylabel(f'clinical patient ({args.patient_col}), ordered by disease axis', fontsize=11)
    # separators between disease states
    prev = None
    for i, st in enumerate(pat['disease_state']):
        if prev is not None and st != prev:
            ax_h.axhline(i - 0.5, color='#333', lw=0.7)
        prev = st
    cb = fig.colorbar(im, ax=ax_h, fraction=0.025, pad=0.02)
    cb.set_label('signature score (z per signature)', fontsize=9)
    cb.ax.tick_params(labelsize=8)
    ax_h.set_title('LTL331 cluster signatures across clinical patients, by disease axis',
                   fontsize=9)
    savefig(fig, args.out, 'clinical_signature_heatmap', args.plot_format)
    print(f'\n[ok] wrote patient_signature_scores.tsv, signature_trend_tests.tsv, '
          f'clinical_signature_heatmap.{{{",".join(args.plot_format)}}} to {args.out}/')


if __name__ == '__main__':
    main()
