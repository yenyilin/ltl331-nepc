#!/usr/bin/env python3
"""
cellrank_bifurcation.py  —  does cluster 17 give rise to BOTH ASCL1+ and
ASCL1- NEPC?  (scVelo + CellRank2; run on the LTL331-only object on the server)

Designed to be able to FALSIFY the bifurcation claim, not just confirm it. It
distinguishes three hypotheses:
  (H1) bifurcation  — cluster 17 feeds both ASCL1+ (cluster 10) and ASCL1- (cl 7)
  (H2) biased       — cluster 17 mostly feeds one fate
  (H3) separate     — the ASCL1- state does not descend from cluster 17

Pipeline:
  0  load LTL331 object (with spliced/unspliced or precomputed velocity)
  1  velocity (optional recompute; dynamical recommended) + velocity graph
     (skipped when --kernel pseudotime/realtime — neither uses scvelo)
  2  kernel: chosen by --kernel
       velocity (default): (1-w)*VelocityKernel   + w*ConnectivityKernel
       pseudotime:         (1-w)*PseudotimeKernel + w*ConnectivityKernel
                            (uses --pseudotime-key, e.g. velocity_pseudotime;
                             use when velocity direction is unreliable but a
                             graph pseudotime tracks real time)
       realtime:           moscot.TemporalProblem -> RealTimeKernel
                            (uses --time-key + entropic OT with --ot-epsilon)
  Run all three into separate --outdir and compare the cluster-17 fate split for
  a kernel-invariance table (the strongest form of the R1 #5 result).
  3  GPCCA macrostates DATA-DRIVEN: do ASCL1+/- emerge as terminals on their own?
  4  DIRECTED fate probabilities with terminals = {ASCL1+ cluster, ASCL1- cluster}
  5  cluster-17 fate-probability distribution + split fraction + paired test
  6  controls: positive (10->ASCL1+, 7->ASCL1-) and negative (PRAD->~0)
  7  robustness: sweep kernel weight + random seed; report stability
  8  lineage drivers (ASCL1 should top the ASCL1+ lineage)
  9  separate-origin check: does ASCL1- receive flow from non-17 routes?

Deps: cellrank>=2.0, scvelo, scanpy/anndata, numpy, pandas, scipy, matplotlib.
Outputs: TSVs + figures (fonttype 42) in --outdir.

Example (velocity kernel):
  python cellrank_bifurcation.py \\
      --h5ad ltl331_velocity.h5ad --outdir cr_bifurcation_vel \\
      --cluster-key seurat_clusters --intermediate 17 \\
      --ascl1-pos 10 --ascl1-neg 7 \\
      --recompute-velocity --velocity-mode dynamical \\
      --conn-weight 0.2 --n-macrostates 12 --seeds 0 1 2

Example (pseudotime kernel — no velocity needed):
  python cellrank_bifurcation.py \\
      --h5ad ltl331_velocity.h5ad --outdir cr_bifurcation_pt \\
      --kernel pseudotime --pseudotime-key velocity_pseudotime \\
      --cluster-key seurat_clusters --intermediate 17 \\
      --ascl1-pos 10 --ascl1-neg 7 --n-macrostates 12 --seeds 0 1 2
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _utils import setup_style, savefig


# --------------------------------------------------------------------------- #
# 0-1  load + velocity
# --------------------------------------------------------------------------- #
def load_and_velocity(args):
    import scanpy as sc
    adata = sc.read_h5ad(args.h5ad)
    print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes')
    if args.cluster_key not in adata.obs:
        raise SystemExit(f'cluster key {args.cluster_key} not in obs')
    # CellRank maps macrostates -> clusters via the .cat accessor, so the cluster key MUST be
    # a pandas categorical (Seurat->h5ad often returns it as object/int -> "Can only use .cat
    # accessor with a 'category' dtype"). Cast to clean string-categories; downstream
    # .astype(str) on a categorical still works.
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str).astype('category')

    if args.kernel == 'realtime':
        # realtime kernel doesn't use scvelo at all; we only need the object
        if args.time_key not in adata.obs.columns:
            raise SystemExit(f'--time-key "{args.time_key}" not in obs '
                             f'(have: {list(adata.obs.columns)})')
        print(f'[realtime] using --time-key={args.time_key}; '
              f'scvelo / velocity_confidence skipped')
        return adata

    if args.kernel == 'pseudotime':
        # pseudotime kernel needs only a precomputed pseudotime + the kNN graph;
        # no RNA-velocity layer required
        if args.pseudotime_key not in adata.obs.columns:
            raise SystemExit(f'--pseudotime-key "{args.pseudotime_key}" not in obs '
                             f'(have: {list(adata.obs.columns)})')
        n_nan = int(adata.obs[args.pseudotime_key].isna().sum())
        if n_nan:
            print(f'[pseudotime] WARNING: {n_nan} NaN values in '
                  f'{args.pseudotime_key}; drop or impute before running '
                  f'(PseudotimeKernel cannot handle NaNs)')
        print(f'[pseudotime] using --pseudotime-key={args.pseudotime_key}; '
              f'scvelo / velocity recompute skipped')
        return adata

    # velocity kernel path
    import scvelo as scv
    if args.recompute_velocity:
        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
            raise SystemExit('spliced/unspliced layers required to recompute '
                             'velocity')
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        if args.velocity_mode == 'dynamical':
            scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
            scv.tl.velocity(adata, mode='dynamical')
        else:
            scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
    else:
        if 'velocity' not in adata.layers:
            raise SystemExit('no velocity layer found; pass --recompute-velocity')

    # reliability gate: velocity confidence in the intermediate cluster
    try:
        scv.tl.velocity_confidence(adata)
        conf = adata.obs.groupby(args.cluster_key)['velocity_confidence'].median()
        conf.to_csv(Path(args.outdir) / 'velocity_confidence_by_cluster.tsv',
                    sep='\t')
        ic = conf.get(str(args.intermediate), float('nan'))
        print(f'[gate] median velocity_confidence in cluster '
              f'{args.intermediate} = {ic:.3f} '
              f'(low => absorption estimates unreliable, report caveat)')
    except Exception as e:
        print(f'[gate] velocity_confidence skipped: {e}')
    return adata


# --------------------------------------------------------------------------- #
# 2  kernel — dispatches on args.kernel
# --------------------------------------------------------------------------- #
def _build_velocity_kernel(adata, conn_weight):
    from cellrank.kernels import VelocityKernel, ConnectivityKernel
    vk = VelocityKernel(adata).compute_transition_matrix()
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    w = 1.0 - conn_weight
    return w * vk + conn_weight * ck


def _build_pseudotime_kernel(adata, pseudotime_key, conn_weight):
    """(1-w)*PseudotimeKernel + w*ConnectivityKernel.

    Directionality comes from a precomputed pseudotime (e.g. velocity_pseudotime);
    no RNA-velocity layer is used. Soft thresholding biases transitions toward
    increasing pseudotime without zeroing out the local neighbourhood, which is
    more robust than the hard scheme near the root.
    """
    from cellrank.kernels import PseudotimeKernel, ConnectivityKernel
    pk = PseudotimeKernel(adata, time_key=pseudotime_key).compute_transition_matrix(
        threshold_scheme='soft')
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    w = 1.0 - conn_weight
    return w * pk + conn_weight * ck


def _build_realtime_kernel(adata, time_key, ot_epsilon):
    """RealTimeKernel via moscot OT between consecutive --time-key timepoints."""
    try:
        from moscot.problems.time import TemporalProblem
        from cellrank.kernels import RealTimeKernel
    except ImportError as e:
        raise SystemExit(
            'moscot / cellrank RealTimeKernel not importable. Install with:\n'
            '    uv pip install "moscot[scrnaseq]"\n'
            f'  ({e})'
        )
    if time_key not in adata.obs.columns:
        raise SystemExit(f'--time-key "{time_key}" not in adata.obs '
                         f'(have: {list(adata.obs.columns)})')
    print(f'[realtime] solving moscot.TemporalProblem on time_key={time_key} '
          f'(epsilon={ot_epsilon}) — this is CPU-bound on Intel Mac and can '
          f'take minutes to hours')
    tp = TemporalProblem(adata)
    tp.prepare(time_key=time_key)
    tp.solve(epsilon=ot_epsilon)
    return RealTimeKernel.from_moscot(tp).compute_transition_matrix()


def build_kernel(adata, args, conn_weight=None):
    """Return a CellRank2 kernel chosen by args.kernel.

      'velocity'   -> (1-w)*VelocityKernel   + w*ConnectivityKernel
      'pseudotime' -> (1-w)*PseudotimeKernel + w*ConnectivityKernel
                      (direction from args.pseudotime_key)
      'realtime'   -> moscot OT between consecutive --time-key timepoints
                      wrapped in RealTimeKernel
    """
    if conn_weight is None:
        conn_weight = args.conn_weight
    if args.kernel == 'velocity':
        return _build_velocity_kernel(adata, conn_weight)
    if args.kernel == 'pseudotime':
        return _build_pseudotime_kernel(adata, args.pseudotime_key, conn_weight)
    if args.kernel == 'realtime':
        return _build_realtime_kernel(adata, args.time_key, args.ot_epsilon)
    raise SystemExit(f'unknown --kernel value: {args.kernel}')


# --------------------------------------------------------------------------- #
# 3  data-driven macrostates (do ASCL1+/- emerge unforced?)
# --------------------------------------------------------------------------- #
def data_driven_macrostates(kernel, adata, args, outdir):
    from cellrank.estimators import GPCCA
    g = GPCCA(kernel)
    g.compute_schur(n_components=max(args.n_macrostates + 4, 20))
    g.compute_macrostates(n_states=args.n_macrostates,
                          cluster_key=args.cluster_key)
    g.predict_terminal_states()

    # which Seurat clusters dominate each terminal macrostate?
    ms = g.macrostates.cat.categories.tolist()
    term = list(g.terminal_states.cat.categories) \
        if g.terminal_states is not None else []
    assign = pd.DataFrame({
        'macrostate': g.macrostates.astype(str),
        'cluster': adata.obs[args.cluster_key].astype(str).values,
    })
    comp = (assign.groupby('macrostate')['cluster']
            .value_counts(normalize=True).rename('frac').reset_index())
    comp.to_csv(outdir / 'macrostate_cluster_composition.tsv', sep='\t',
                index=False)
    with open(outdir / 'macrostates_summary.txt', 'w') as fh:
        fh.write(f'n_macrostates={args.n_macrostates}\n')
        fh.write(f'all macrostates: {ms}\n')
        fh.write(f'predicted terminal states: {term}\n')
    # did ASCL1+ (pos) and ASCL1- (neg) each dominate a terminal macrostate?
    pos = [str(c) for c in args.ascl1_pos]
    neg = [str(c) for c in args.ascl1_neg]
    dom = comp.sort_values('frac', ascending=False).groupby('macrostate').first()
    pos_term = dom[dom['cluster'].isin(pos)].index.tolist()
    neg_term = dom[dom['cluster'].isin(neg)].index.tolist()
    print(f'[macrostates] terminal states (data-driven): {term}')
    print(f'[macrostates] macrostate(s) dominated by ASCL1+ cluster {pos}: {pos_term}')
    print(f'[macrostates] macrostate(s) dominated by ASCL1- cluster {neg}: {neg_term}')
    print('[macrostates] If BOTH appear as terminals unforced -> supports H1 '
          'before any forcing.')
    return g


# --------------------------------------------------------------------------- #
# 4-5  directed fate probabilities + cluster-17 analysis
# --------------------------------------------------------------------------- #
def directed_fates(kernel, adata, args, outdir, seed='main', random_state=None):
    """`seed` is the output-filename tag; `random_state`, if set, re-seeds
    np.random before the GPCCA calls so pygpcca's PCCA+ K-means is reproducible.
    The Schur decomposition itself is deterministic given the transition matrix,
    so observable variability across random_state values is small but non-zero.
    """
    if random_state is not None:
        np.random.seed(int(random_state))
    from cellrank.estimators import GPCCA
    g = GPCCA(kernel)
    g.compute_schur(n_components=max(args.n_macrostates + 4, 20))
    g.compute_macrostates(n_states=args.n_macrostates,
                          cluster_key=args.cluster_key)

    # set terminal states to the two NE fates (each may be a GROUP of clusters)
    pos = [str(c) for c in args.ascl1_pos]
    neg = [str(c) for c in args.ascl1_neg]
    try:
        g.set_terminal_states({'ASCL1pos': _cells_of(adata, args, pos),
                               'ASCL1neg': _cells_of(adata, args, neg)})
    except Exception:
        # fall back to macrostate-name based selection
        g.set_terminal_states([m for m in g.macrostates.cat.categories
                               if any(p in m for p in pos) or any(n in m for n in neg)])
    g.compute_fate_probabilities()
    fp = g.fate_probabilities
    fp_df = pd.DataFrame(fp.X, index=adata.obs_names, columns=fp.names)
    fp_df[args.cluster_key] = adata.obs[args.cluster_key].astype(str).values
    fp_df.to_csv(outdir / f'fate_probabilities_seed{seed}.tsv', sep='\t')
    return g, fp_df


def _cells_of(adata, args, clusters):
    """Cells belonging to any of `clusters` (a str or a list of str)."""
    if isinstance(clusters, str):
        clusters = [clusters]
    clset = {str(c) for c in clusters}
    return adata.obs_names[adata.obs[args.cluster_key].astype(str).isin(clset)]


def analyze_intermediate(fp_df, args, outdir, plot_formats):
    pos_col = [c for c in fp_df.columns if 'pos' in c.lower()]
    neg_col = [c for c in fp_df.columns if 'neg' in c.lower()]
    if not pos_col or not neg_col:
        print(f'[intermediate] could not find pos/neg fate columns in '
              f'{list(fp_df.columns)}')
        return
    pcol, ncol = pos_col[0], neg_col[0]
    inter = fp_df[fp_df[args.cluster_key] == str(args.intermediate)]
    if inter.empty:
        print(f'[intermediate] no cells in cluster {args.intermediate}')
        return
    p_pos = inter[pcol].values
    p_neg = inter[ncol].values
    diff = p_pos - p_neg
    frac_pos = float((p_pos > p_neg).mean())
    # paired Wilcoxon (descriptive within one PDX; cells not independent)
    try:
        w, pval = stats.wilcoxon(p_pos, p_neg)
    except Exception:
        w, pval = float('nan'), float('nan')
    # baseline root split (PRAD/root clusters): absorption probs sum to 1, so the
    # meaningful quantity is the intermediate's ENRICHMENT over this baseline,
    # not its value relative to 0
    prad = [str(c) for c in args.prad_clusters]
    base = fp_df[fp_df[args.cluster_key].isin(prad)]
    base_neg = float(base[ncol].mean()) if not base.empty else float('nan')
    summary = {
        'n_cluster17_cells': int(len(inter)),
        'mean_P_ASCL1pos': float(np.mean(p_pos)),
        'mean_P_ASCL1neg': float(np.mean(p_neg)),
        'median_diff_pos_minus_neg': float(np.median(diff)),
        'frac_biased_to_ASCL1pos': frac_pos,
        'frac_biased_to_ASCL1neg': 1 - frac_pos,
        'baseline_root_P_ASCL1neg': base_neg,
        'enrichment_P_ASCL1neg_over_baseline': float(np.mean(p_neg)) - base_neg,
        'wilcoxon_W': float(w), 'wilcoxon_p': float(pval),
    }
    pd.Series(summary).to_csv(outdir / 'intermediate_fate_summary.tsv', sep='\t')
    print('[intermediate] cluster-17 fate summary:')
    for k, v in summary.items():
        print(f'    {k}: {v}')
    print('    INTERPRET: mean P_pos ~= P_neg and both substantial -> H1 '
          '(bifurcation). Strong skew -> H2 (biased). Near-zero to one fate '
          '-> that fate does not derive from cluster 17 (toward H3). Judge bias '
          'as ENRICHMENT over baseline_root_P_ASCL1neg (probs sum to 1; the root '
          'is not 0), e.g. P_neg 0.87 vs baseline 0.40 = strong ASCL1- gateway.')

    plt = setup_style()
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.violinplot([p_pos, p_neg], showmeans=True)
    ax.set_xticks([1, 2]); ax.set_xticklabels(['ASCL1+', 'ASCL1-'])
    ax.set_ylabel('absorption probability'); ax.set_ylim(0, 1)
    ax.set_title(f'Cluster {args.intermediate} fate (n={len(inter)})')
    savefig(fig, outdir, 'intermediate_fate_violin', plot_formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# 10  per-cluster fates (H3a: separate route?) — does any non-intermediate
#     cluster also have substantial P(ASCL1-) without going through 17?
# --------------------------------------------------------------------------- #
def upstream_fates(fp_df, args, outdir, plot_formats):
    """Per-cluster mean fate probabilities, sorted by P(ASCL1-).

    Tests H3a (separate route to ASCL1-): if cluster 17 has the highest
    P(ASCL1-) among non-terminal clusters, it is the primary route. If some
    non-17 cluster has comparable/higher P(ASCL1-), a bypass route exists —
    the bifurcation is multi-channel rather than exclusively via 17.
    """
    pcol = [c for c in fp_df.columns if 'pos' in c.lower()][0]
    ncol = [c for c in fp_df.columns if 'neg' in c.lower()][0]
    per = (fp_df.groupby(args.cluster_key)[[pcol, ncol]].mean()
                .rename(columns={pcol: 'mean_P_pos', ncol: 'mean_P_neg'}))
    thr = args.upstream_threshold
    per['both_substantial'] = (per['mean_P_pos'] > thr) & (per['mean_P_neg'] > thr)
    def _bias(r):
        if r['mean_P_pos'] - r['mean_P_neg'] > 0.2: return 'ASCL1+ biased'
        if r['mean_P_neg'] - r['mean_P_pos'] > 0.2: return 'ASCL1- biased'
        return 'mixed'
    per['bias'] = per.apply(_bias, axis=1)
    per = per.sort_values('mean_P_neg', ascending=False)
    per.to_csv(outdir / 'step10_per_cluster_fates.tsv', sep='\t')

    inter = str(args.intermediate)
    print('[step10] per-cluster mean fate probabilities (sorted by P(ASCL1-)):')
    print(per.round(3).to_string())
    if inter in per.index:
        rank = int(per['mean_P_neg'].rank(ascending=False).loc[inter])
        print(f'[step10] cluster {inter} ranks {rank}/{len(per)} by mean P(ASCL1-)')
        higher = per[per['mean_P_neg'] > per.loc[inter, 'mean_P_neg']].index.tolist()
        neg_set = {str(c) for c in args.ascl1_neg}
        higher = [c for c in higher if c not in neg_set]  # drop the terminal clusters
        if higher:
            print(f'[step10] non-{inter} clusters with higher mean P(ASCL1-) '
                  f'than {inter}: {higher} -- possible bypass route(s) '
                  f'(H3a partially supported; the trajectory is multi-channel)')
        else:
            print(f'[step10] no non-terminal cluster outranks {inter} for P(ASCL1-) '
                  f'-- cluster {inter} is the primary route to ASCL1- (H3a not supported)')

    plt = setup_style()
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(per['mean_P_pos'], per['mean_P_neg'], s=30,
               edgecolor='k', linewidth=0.4, alpha=0.85)
    for cl, row in per.iterrows():
        is_inter = (cl == inter)
        ax.annotate(str(cl),
                    (row['mean_P_pos'], row['mean_P_neg']),
                    fontsize=6,
                    color='red' if is_inter else 'black',
                    fontweight='bold' if is_inter else 'normal',
                    xytext=(3, 3), textcoords='offset points')
    ax.axhline(thr, color='grey', lw=0.4, ls='--')
    ax.axvline(thr, color='grey', lw=0.4, ls='--')
    ax.set_xlabel('mean P(ASCL1+)'); ax.set_ylabel('mean P(ASCL1-)')
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_title(f'Per-cluster mean fate probabilities (intermediate={inter} red)')
    savefig(fig, outdir, 'step10_cluster_fate_scatter', plot_formats)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# 6  controls
# --------------------------------------------------------------------------- #
def controls(fp_df, args, outdir):
    pcol = [c for c in fp_df.columns if 'pos' in c.lower()][0]
    ncol = [c for c in fp_df.columns if 'neg' in c.lower()][0]
    rows = []
    # POSITIVE controls: each terminal cluster should self-absorb to its own
    # fate (~1). This is the only true control here.
    checks = {}
    for cl in [str(c) for c in args.ascl1_pos]:
        checks[f'POSctrl_cluster{cl}_to_ASCL1pos'] = (cl, pcol, '~1 (terminal self-absorption)')
    for cl in [str(c) for c in args.ascl1_neg]:
        checks[f'POSctrl_cluster{cl}_to_ASCL1neg'] = (cl, ncol, '~1 (terminal self-absorption)')
    # BASELINE (NOT a "~0" control): PRAD/root clusters are ON the trajectory, so
    # their fate probs are NOT expected to be 0 -- absorption probabilities sum to
    # 1 per cell, so the root necessarily distributes across the terminals. These
    # rows are the baseline root split, the reference against which the
    # intermediate cluster's enrichment is judged.
    for cl in [str(c) for c in args.prad_clusters]:
        checks[f'BASELINE_root{cl}_ASCL1pos'] = (cl, pcol, 'root split (sums to 1; reference)')
        checks[f'BASELINE_root{cl}_ASCL1neg'] = (cl, ncol, 'root split (sums to 1; reference)')
    for name, (cl, col, expect) in checks.items():
        sub = fp_df[fp_df[args.cluster_key] == cl]
        val = float(sub[col].mean()) if not sub.empty else float('nan')
        rows.append({'check': name, 'cluster': cl, 'fate': col,
                     'mean_prob': val, 'expected': expect})
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'controls.tsv', sep='\t', index=False)
    print('[controls] terminals should self-absorb (~1); PRAD/root rows are the '
          'BASELINE split (sum to 1), NOT expected ~0:')
    print(df.to_string(index=False))

    # baseline ASCL1- fraction (mean over PRAD/root clusters): the reference the
    # intermediate cluster should be compared against (enrichment, not vs 0)
    prad = [str(c) for c in args.prad_clusters]
    base = fp_df[fp_df[args.cluster_key].isin(prad)]
    if not base.empty:
        print(f'[controls] baseline root P(ASCL1-) = {float(base[ncol].mean()):.3f} '
              f'(mean over PRAD {prad}). Judge the intermediate\'s P(ASCL1-) as '
              f'ENRICHMENT over this baseline, not relative to 0.')


# --------------------------------------------------------------------------- #
# 7  robustness across seeds / kernel weights
# --------------------------------------------------------------------------- #
def robustness(adata, args, outdir):
    if args.kernel == 'realtime':
        print(f'[robustness] skipped: the conn-weight sweep mixes the primary '
              f'kernel with ConnectivityKernel, which the moscot RealTimeKernel '
              f'does not use. For realtime sensitivity, rerun with a different '
              f'--ot-epsilon.')
        return
    # velocity and pseudotime both mix in ConnectivityKernel via conn_weight,
    # so the sweep is meaningful for both
    rows = []
    for w in args.conn_weights:
        for seed in args.seeds:
            # Defensive copy: GPCCA's set_terminal_states / compute_macrostates
            # write to adata.obs/uns. Reusing the same adata across sweep
            # iterations leaks residual state into the next solve.
            adata_iter = adata.copy()
            kernel = build_kernel(adata_iter, args, conn_weight=w)
            try:
                _, fp_df = directed_fates(kernel, adata_iter, args, outdir,
                                          seed=f'w{w}_s{seed}',
                                          random_state=seed)
                pcol = [c for c in fp_df.columns if 'pos' in c.lower()][0]
                ncol = [c for c in fp_df.columns if 'neg' in c.lower()][0]
                inter = fp_df[fp_df[args.cluster_key] == str(args.intermediate)]
                rows.append({'conn_weight': w, 'seed': seed,
                             'mean_P_pos': float(inter[pcol].mean()),
                             'mean_P_neg': float(inter[ncol].mean())})
            except Exception as e:
                print(f'[robustness] w={w} seed={seed} failed: {e}')
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'robustness_sweep.tsv', sep='\t', index=False)
    if not df.empty:
        print('[robustness] cluster-17 mean fate across params:')
        print(df.to_string(index=False))
        print(f'    P_pos range: {df.mean_P_pos.min():.3f}-{df.mean_P_pos.max():.3f}; '
              f'P_neg range: {df.mean_P_neg.min():.3f}-{df.mean_P_neg.max():.3f} '
              f'(tight ranges => stable bifurcation call)')


# --------------------------------------------------------------------------- #
# SAVE — persist the canonical run (terminal states + fate probabilities +
# parameters) into an h5ad + a clean TSV. Call right after the canonical run
# and BEFORE the robustness sweep (which mutates adata.obsp).
# --------------------------------------------------------------------------- #
def save_canonical(g, fp_df, adata, args, outdir):
    """g.write_to_adata() + uns params + fate-probs TSV.

    The h5ad is fully reloadable:
        adata = anndata.read_h5ad(<path>)
        g = cellrank.estimators.GPCCA.read_from_adata(adata)
    """
    import cellrank as cr_mod
    try:
        import scvelo as scv_mod
        scv_ver = getattr(scv_mod, '__version__', 'unknown')
    except Exception:
        scv_ver = 'not-imported'

    g.write_to_adata()    # writes macrostates, terminal_states, fate_probabilities

    if args.kernel == 'velocity':
        kernel_descr = (f'{1.0 - args.conn_weight:.2f}*VelocityKernel + '
                        f'{args.conn_weight:.2f}*ConnectivityKernel')
    elif args.kernel == 'pseudotime':
        kernel_descr = (f'{1.0 - args.conn_weight:.2f}*PseudotimeKernel + '
                        f'{args.conn_weight:.2f}*ConnectivityKernel '
                        f'(pseudotime_key={args.pseudotime_key})')
    else:
        kernel_descr = (f'RealTimeKernel via moscot.TemporalProblem '
                        f'(time_key={args.time_key}, epsilon={args.ot_epsilon})')

    adata.uns['cellrank_params'] = {
        'cellrank_version': getattr(cr_mod, '__version__', 'unknown'),
        'scvelo_version': scv_ver,
        'kernel_type': args.kernel,
        'kernel': kernel_descr,
        'n_macrostates': int(args.n_macrostates),
        'cluster_key': args.cluster_key,
        'intermediate_cluster': str(args.intermediate),
        'ascl1_pos_cluster': ','.join(str(c) for c in args.ascl1_pos),
        'ascl1_neg_cluster': ','.join(str(c) for c in args.ascl1_neg),
        'prad_negative_control_clusters': [str(c) for c in args.prad_clusters],
        'velocity_mode': ((args.velocity_mode if args.recompute_velocity
                           else 'as-provided') if args.kernel == 'velocity'
                          else f'n/a ({args.kernel} kernel)'),
        'pseudotime_key': (args.pseudotime_key if args.kernel == 'pseudotime'
                           else 'n/a'),
        'note': ('Canonical bifurcation analysis. Reload with '
                 'cellrank.estimators.GPCCA.read_from_adata(adata).'),
    }

    suffix = args.save_suffix or (
        f"pos{'-'.join(str(c) for c in args.ascl1_pos)}"
        f"_neg{'-'.join(str(c) for c in args.ascl1_neg)}")
    h5_path = Path(outdir) / f'cellrank_canonical_{suffix}.h5ad'
    print(f'[save] writing {h5_path} (gzip)')
    adata.write_h5ad(h5_path, compression='gzip')

    tsv = fp_df.drop(columns=[args.cluster_key], errors='ignore').copy()
    tsv_path = Path(outdir) / f'fate_probabilities_canonical_{suffix}.tsv'
    tsv.to_csv(tsv_path, sep='\t')

    print(f'[save] wrote {h5_path}')
    print(f'[save] wrote {tsv_path}')
    print(f'[save] reload: adata = anndata.read_h5ad("{h5_path.name}"); '
          f'g = cellrank.estimators.GPCCA.read_from_adata(adata)')


# --------------------------------------------------------------------------- #
# 8  lineage drivers (confirm terminals are really ASCL1+/-)
# --------------------------------------------------------------------------- #
def lineage_drivers(g, outdir):
    try:
        drivers = g.compute_lineage_drivers()
        drivers.to_csv(outdir / 'lineage_drivers.tsv', sep='\t')
        print('[drivers] wrote lineage_drivers.tsv — check ASCL1 ranks top for '
              'the ASCL1+ lineage and differs for ASCL1-')
    except Exception as e:
        print(f'[drivers] skipped: {e}')


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--outdir', default='cr_bifurcation')
    p.add_argument('--cluster-key', default='clusters')
    p.add_argument('--intermediate', default='17')
    p.add_argument('--ascl1-pos', nargs='+', default=['10'],
                   help='ASCL1+ terminal cluster(s); space-separated for a group, '
                        'e.g. --ascl1-pos 10 14')
    p.add_argument('--ascl1-neg', nargs='+', default=['7'],
                   help='ASCL1- terminal cluster(s); space-separated for a group, '
                        'e.g. --ascl1-neg 1 3 9')
    p.add_argument('--prad-clusters', nargs='+', default=['11', '15'],
                   help='AR-high PRAD clusters for the negative control')
    p.add_argument('--kernel', default='velocity',
                   choices=['velocity', 'pseudotime', 'realtime'],
                   help='CellRank2 kernel: velocity (VelocityKernel + '
                        'ConnectivityKernel, canonical scvelo recipe); '
                        'pseudotime (PseudotimeKernel on --pseudotime-key + '
                        'ConnectivityKernel — use when velocity direction is '
                        'unreliable but a graph pseudotime tracks real time); '
                        'or realtime (RealTimeKernel via moscot.TemporalProblem '
                        'using --time-key + entropic OT). Run each into a '
                        'separate --outdir for the kernel-invariance check.')
    p.add_argument('--time-key', default='timepoint',
                   help='obs column with experimental timepoints (used by '
                        '--kernel realtime); e.g. preCx, 2wk, 4wk, ...')
    p.add_argument('--pseudotime-key', default='velocity_pseudotime',
                   help='obs column with a precomputed pseudotime (used by '
                        '--kernel pseudotime); e.g. velocity_pseudotime, '
                        'dpt_pseudotime. Must be NaN-free.')
    p.add_argument('--ot-epsilon', type=float, default=1e-2,
                   help='moscot OT entropic regularisation '
                        '(--kernel realtime only)')
    p.add_argument('--recompute-velocity', action='store_true')
    p.add_argument('--velocity-mode', default='dynamical',
                   choices=['dynamical', 'stochastic'])
    p.add_argument('--conn-weight', type=float, default=0.2,
                   help='weight of ConnectivityKernel in the main run')
    p.add_argument('--conn-weights', nargs='+', type=float,
                   default=[0.1, 0.2, 0.3], help='sweep for robustness')
    p.add_argument('--n-macrostates', type=int, default=12)
    p.add_argument('--seeds', nargs='+', type=int, default=[0, 1, 2])
    p.add_argument('--n-jobs', type=int, default=8)
    p.add_argument('--upstream-threshold', type=float, default=0.2,
                   help='P threshold above which a cluster is treated as having '
                        'substantial access to a fate (Step 10, H3a check)')
    p.add_argument('--save-adata', action='store_true',
                   help='persist the canonical run (terminals + fate probs + '
                        'parameters) to a gzipped h5ad and a TSV')
    p.add_argument('--save-suffix', default=None,
                   help='filename suffix (default: pos<X>_neg<Y> from terminals)')
    p.add_argument('--plot-format', nargs='+', default=['pdf', 'png'],
                   choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    args.outdir = str(outdir)

    adata = load_and_velocity(args)
    kernel = build_kernel(adata, args)

    print('\n=== STEP 3: data-driven macrostates (unforced) ===')
    g_dd = data_driven_macrostates(kernel, adata, args, outdir)

    print('\n=== STEP 4-5: directed fate probabilities + cluster-17 analysis ===')
    g_dir, fp_df = directed_fates(kernel, adata, args, outdir, seed='main',
                                  random_state=args.seeds[0])
    analyze_intermediate(fp_df, args, outdir, args.plot_format)

    print('\n=== STEP 6: controls ===')
    controls(fp_df, args, outdir)

    print('\n=== STEP 10: per-cluster fates / bypass-route check (H3a) ===')
    upstream_fates(fp_df, args, outdir, args.plot_format)

    print('\n=== STEP 8: lineage drivers ===')
    lineage_drivers(g_dir, outdir)

    if args.save_adata:
        print('\n=== SAVE: canonical h5ad + TSV ===')
        save_canonical(g_dir, fp_df, adata, args, outdir)

    print('\n=== STEP 7: robustness sweep ===')
    robustness(adata, args, outdir)

    print('\nDone. Review intermediate_fate_summary.tsv + controls.tsv + '
          'robustness_sweep.tsv together to decide H1/H2/H3.')


if __name__ == '__main__':
    main()
