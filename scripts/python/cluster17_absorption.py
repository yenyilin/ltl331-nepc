#!/usr/bin/env python3
"""
cluster17_absorption.py — absorption probabilities recomputed on THIS object's
actual CellRank terminal macrostates (reproducible; companion to the transition-flux).

The object's terminal states (obs['term_states_fwd']) are {0,1,4_1,4_2,6,7,9,10} — so
the ASCL1+ state (cluster 10) IS a recovered terminal here, and the NE compartment
comprises several terminal macrostates (1,7,9,10). We recompute per-cell fate (absorption)
probabilities with the SAME kernel used elsewhere
(PseudotimeKernel(velocity_pseudotime)*0.8 + ConnectivityKernel*0.2), set the terminal
states to the object's own term_states_fwd, and summarize the fate of cluster-17 cells:

  - P(NE)  = sum of fate mass over NE terminals      {1,7,9,10}
  - P(nonNE) = sum over non-NE terminals             {0,4_1,4_2,6}
  - within the NE compartment: P(ASCL1+)= {10} vs P(ASCL1-)= {1,7,9}

Reports median/mean over the cluster-17 cells + a paired Wilcoxon signed-rank test
(P_NE > P_nonNE). These numbers are reproducible from the shared object.

Example:
  python cluster17_absorption.py --h5ad ltl331_velocity_cellrank.h5ad \\
      --cluster-key seurat_clusters --source 17 \\
      --term-key term_states_fwd --pt-key velocity_pseudotime --conn-weight 0.2 \\
      --ne 1 7 9 10 --non-ne 0 4_1 4_2 6 --ascl1-pos 10 --ascl1-neg 1 7 9 \\
      --outdir data/cluster17_absorption
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def grp(fp_df, names):
    """sum fate-probability columns whose terminal name is in `names` (string match)."""
    cols = [c for c in fp_df.columns if str(c) in {str(n) for n in names}]
    missing = {str(n) for n in names} - {str(c) for c in fp_df.columns}
    if missing:
        print(f'[warn] terminals not found in fate_probabilities and skipped: {sorted(missing)}')
    return fp_df[cols].sum(axis=1), cols


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--cluster-key', default='seurat_clusters')
    ap.add_argument('--source', default='17')
    ap.add_argument('--term-key', default='term_states_fwd')
    ap.add_argument('--pt-key', default='velocity_pseudotime')
    ap.add_argument('--conn-weight', type=float, default=0.2)
    ap.add_argument('--ne', nargs='+', default=['1', '7', '9', '10'])
    ap.add_argument('--non-ne', nargs='+', default=['0', '4_1', '4_2', '6'])
    ap.add_argument('--ascl1-pos', nargs='+', default=['10'])
    ap.add_argument('--ascl1-neg', nargs='+', default=['1', '7', '9'])
    ap.add_argument('--outdir', default='cluster17_absorption')
    args = ap.parse_args()

    import anndata as ad
    import cellrank as cr
    from scipy.stats import wilcoxon

    adata = ad.read_h5ad(args.h5ad)
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    # g.set_terminal_states reads .cat.categories on the terminal column, so it MUST be
    # categorical (Seurat->h5ad often returns it object/str -> "Can only use .cat accessor
    # with a 'category' dtype"). astype('category') preserves NaN for non-terminal cells.
    adata.obs[args.term_key] = adata.obs[args.term_key].astype('category')
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print('terminal states in obs[%r]:' % args.term_key)
    print(adata.obs[args.term_key].value_counts().to_string())

    # same kernel as cluster17_transition_flux.py --recompute
    pk = cr.kernels.PseudotimeKernel(adata, time_key=args.pt_key).compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    comb = (1 - args.conn_weight) * pk + args.conn_weight * ck
    if getattr(comb, 'transition_matrix', None) is None:
        comb.compute_transition_matrix()

    g = cr.estimators.GPCCA(comb)
    g.set_terminal_states(adata.obs[args.term_key])         # use the object's own terminals
    g.compute_fate_probabilities()
    fp = g.fate_probabilities
    fp_df = pd.DataFrame(np.asarray(fp.X), columns=[str(n) for n in fp.names],
                         index=adata.obs_names)
    fp_df.to_csv(outdir / 'fate_probabilities_allcells.tsv', sep='\t')

    src = adata.obs[args.cluster_key].values == str(args.source)
    n_src = int(src.sum())
    sub = fp_df.loc[src]

    p_ne, ne_cols = grp(sub, args.ne)
    p_non, non_cols = grp(sub, args.non_ne)
    p_pos, _ = grp(sub, args.ascl1_pos)
    p_neg, _ = grp(sub, args.ascl1_neg)
    # within-NE normalized ASCL1 split
    ne_tot = (p_pos + p_neg).replace(0, np.nan)
    frac_pos = (p_pos / ne_tot)

    print(f'\n=== cluster {args.source} fate (n={n_src}); NE terminals={ne_cols}; '
          f'non-NE={non_cols} ===')
    def s(x): return f'median={np.nanmedian(x):.4f}  mean={np.nanmean(x):.4f}'
    print(f'  P(NE)            : {s(p_ne)}')
    print(f'  P(non-NE)        : {s(p_non)}')
    print(f'  P(ASCL1+ {args.ascl1_pos}) : {s(p_pos)}')
    print(f'  P(ASCL1- {args.ascl1_neg}) : {s(p_neg)}')
    print(f'  within-NE fraction ASCL1+ : {s(frac_pos)}')
    # NE-vs-non-NE: kept only as a sanity check. NOTE it is weak/attackable — the non-NE
    # terminals (0/4/6) are UPSTREAM clusters labelled terminal by GPCCA, so P(non-NE) is
    # mostly backward-diffusion leakage from the symmetric ConnectivityKernel, not biology.
    p_w_ne_pval = p_w_asc_pval = float('nan')
    try:
        stat, p_w_ne_pval = wilcoxon(p_ne.values, p_non.values, alternative='greater')
        print(f'  [sanity, weak] Wilcoxon P(NE) > P(non-NE): stat={stat:.1f}, P={p_w_ne_pval:.3e}')
    except ValueError as e:
        print(f'  [wilcoxon NE-vs-nonNE skipped] {e}')
    # the test that actually answers the bifurcation: within-NE, forward terminals only
    try:
        stat2, p_w_asc_pval = wilcoxon(p_neg.values, p_pos.values, alternative='greater')
        print(f'  *** Wilcoxon P(ASCL1-) > P(ASCL1+): stat={stat2:.1f}, P={p_w_asc_pval:.3e} '
              f'(the bifurcation test) ***')
    except ValueError as e:
        print(f'  [wilcoxon ASCL1-vs+ skipped] {e}')

    pd.DataFrame({'source': [args.source], 'n_cells': [n_src],
                  'median_P_NE': [float(np.nanmedian(p_ne))],
                  'median_P_nonNE': [float(np.nanmedian(p_non))],
                  'median_P_ASCL1pos': [float(np.nanmedian(p_pos))],
                  'median_P_ASCL1neg': [float(np.nanmedian(p_neg))],
                  'median_withinNE_frac_ASCL1pos': [float(np.nanmedian(frac_pos))],
                  'wilcoxon_P_ASCL1neg_gt_pos': [p_w_asc_pval],
                  'wilcoxon_P_NE_gt_nonNE': [p_w_ne_pval]}
                 ).to_csv(outdir / f'cl{args.source}_fate_summary.tsv', sep='\t', index=False)
    print(f'\n[ok] wrote fate summary + all-cell fate probabilities to {outdir}/')


if __name__ == '__main__':
    main()
