#!/usr/bin/env python3
"""
patient_level_stats.py  —  Patient-level statistics + label-transfer confidence
for the LTL331 + Gao + Li integration, making the clinical co-localization claims
rigorous (R1 #3, R1 #2, R2 #7a, R1 #5 clinical angle).

CORE PRINCIPLE: the unit of replication is the PATIENT (~10), not the cell.
Every test aggregates within patient first, then compares across patients.

Steps (select with --steps, default all):
  1  label transfer + confidence + cross-validation (+ optional bidirectional)
  2  per-patient proportions + k/N-patients presence  (the workhorse plot)
  3  patient-aware differential abundance
  4  pseudo-bulk signature scoring + random-gene-set null
  5  per-patient MNN fraction (co-localization) + optional LISI
  6  disease-axis stratified heatmap + Spearman/Kendall trend tests
     (OPTIONAL — needs --patient-metadata; otherwise skipped cleanly)

Hard deps:  anndata, numpy, pandas, scipy, scikit-learn, matplotlib, seaborn,
            statsmodels
Optional (auto-used if importable, else self-contained fallback / skip):
            celltypist (Step 1), decoupler (Step 4), pertpy/scCODA (Step 3),
            scib (Step 5)

Example:
  python patient_level_stats.py \\
      --h5ad integrated.h5ad --outdir patient_stats \\
      --latent-key X_harmony_dataset \\
      --group-col dataset --ref-label LTL331 \\
      --patient-col patient --ref-cluster-col original_cluster \\
      --target-cluster 17 \\
      --signatures ltl331_cluster_markers.tsv \\
      --da-group-col patient_phenotype \\
      --steps 1 2 3 4 5 --plot-format pdf png

NOTE: edit the CLI defaults to match your obs keys. Steps degrade gracefully:
a step that needs a column you didn't supply prints a skip message and continues.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist  # noqa: F401  (kept for extension)

from _utils import setup_style, savefig, to_dense_1d  # noqa: F401


# --------------------------------------------------------------------------- #
# utilities
# --------------------------------------------------------------------------- #
def have(mod):
    import importlib
    try:
        importlib.import_module(mod)
        return True
    except Exception:
        return False


def load_adata(path):
    try:
        import anndata as ad
    except Exception as e:
        raise SystemExit(f'anndata required to read h5ad ({e}). '
                         f'Activate the env that has scanpy/anndata.')
    return ad.read_h5ad(path)


# --------------------------------------------------------------------------- #
# Step 1 — label transfer + confidence
# --------------------------------------------------------------------------- #
def step1_label_transfer(adata, args, outdir):
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import classification_report

    if args.latent_key not in adata.obsm:
        print(f'[step1] SKIP: latent key {args.latent_key} not in obsm '
              f'({list(adata.obsm)})')
        return
    obs = adata.obs
    Z = np.asarray(adata.obsm[args.latent_key])

    ref_mask = (obs[args.group_col].astype(str) == args.ref_label).values
    qry_mask = ~ref_mask
    if ref_mask.sum() == 0:
        print(f'[step1] SKIP: no reference cells (group=={args.ref_label})')
        return

    y_ref = obs.loc[ref_mask, args.ref_cluster_col].astype(str).values
    Z_ref, Z_qry = Z[ref_mask], Z[qry_mask]

    # --- reference validation: stratified 5-fold CV, per-class recall ---
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    knn = KNeighborsClassifier(n_neighbors=args.knn_k, weights='distance')
    y_cv = cross_val_predict(knn, Z_ref, y_ref, cv=skf)
    rep = classification_report(y_ref, y_cv, output_dict=True, zero_division=0)
    pd.DataFrame(rep).T.to_csv(outdir / 'step1_cv_report.tsv', sep='\t')
    tgt = str(args.target_cluster)
    if tgt in rep:
        print(f'[step1] CV recall for cluster {tgt}: '
              f'{rep[tgt]["recall"]:.3f} (n={int(rep[tgt]["support"])}) '
              f'-- if low, rely on Steps 4-5 for this cluster')

    # --- transfer to clinical cells (primary: distance-weighted kNN) ---
    knn.fit(Z_ref, y_ref)
    proba = knn.predict_proba(Z_qry)
    classes = knn.classes_
    pred = classes[proba.argmax(1)]
    conf = proba.max(1)

    qry = obs.loc[qry_mask, [args.group_col]].copy()
    if args.patient_col in obs.columns:
        qry[args.patient_col] = obs.loc[qry_mask, args.patient_col].values
    qry['pred_cluster'] = pred
    qry['confidence'] = conf
    qry.to_csv(outdir / 'step1_clinical_assignments.tsv', sep='\t')

    # --- threshold sensitivity table ---
    rows = []
    for thr in (0.3, 0.5, 0.7):
        keep = conf >= thr
        n_tgt = int(((pred == tgt) & keep).sum())
        rows.append({'threshold': thr, 'n_assigned_total': int(keep.sum()),
                     f'n_cluster_{tgt}': n_tgt})
    pd.DataFrame(rows).to_csv(outdir / 'step1_threshold_sensitivity.tsv',
                             sep='\t', index=False)
    print('[step1] wrote CV report, clinical assignments, threshold sensitivity')

    # --- optional celltypist concordance ---
    if args.use_celltypist and have('celltypist'):
        try:
            _celltypist_concordance(adata, ref_mask, qry_mask, y_ref, outdir, args)
        except Exception as e:
            print(f'[step1] celltypist concordance skipped: {e}')
    elif args.use_celltypist:
        print('[step1] celltypist requested but not importable; skipping')


def _celltypist_concordance(adata, ref_mask, qry_mask, y_ref, outdir, args):
    import celltypist
    from celltypist import models
    ref = adata[ref_mask].copy()
    ref.obs['_label'] = y_ref
    model = celltypist.train(ref, labels='_label', n_jobs=4,
                             use_SGD=True, feature_selection=True)
    pred = celltypist.annotate(adata[qry_mask].copy(), model=model,
                               majority_voting=False)
    pred.predicted_labels.to_csv(outdir / 'step1_celltypist_pred.tsv', sep='\t')
    print('[step1] wrote celltypist predictions (concordance vs kNN)')


# --------------------------------------------------------------------------- #
# Step 2 — per-patient proportions + presence
# --------------------------------------------------------------------------- #
def step2_proportions(adata, args, outdir, plot_formats):
    assign_path = outdir / 'step1_clinical_assignments.tsv'
    if not assign_path.exists():
        print('[step2] SKIP: run step 1 first (needs clinical assignments)')
        return
    df = pd.read_csv(assign_path, sep='\t', index_col=0)
    if args.patient_col not in df.columns:
        print(f'[step2] SKIP: patient column {args.patient_col} not available')
        return

    # pred_cluster round-trips through CSV; pandas infers int when every label
    # is a digit, leaving `'17' in ct.columns` (with str tgt) silently False.
    # Force string here so target-cluster lookups below behave consistently.
    df['pred_cluster'] = df['pred_cluster'].astype(str)

    df = df[df['confidence'] >= args.conf_threshold]
    # per-patient proportion of each predicted cluster
    ct = pd.crosstab(df[args.patient_col], df['pred_cluster'])
    prop = ct.div(ct.sum(1), axis=0)
    prop.to_csv(outdir / 'step2_per_patient_proportions.tsv', sep='\t')

    # k/N presence for the target cluster
    tgt = str(args.target_cluster)
    n_pat = df[args.patient_col].nunique()
    if tgt in ct.columns:
        per_pat = ct[tgt]
        k = int((per_pat > 0).sum())
        print(f'[step2] cluster {tgt} present in {k}/{n_pat} patients '
              f'(counts: {dict(per_pat)})')
        per_pat.to_csv(outdir / f'step2_cluster{tgt}_per_patient.tsv', sep='\t')
    else:
        print(f'[step2] cluster {tgt} not assigned to any clinical cell '
              f'at confidence >= {args.conf_threshold}')

    # plot: stacked per-patient composition
    plt = setup_style()
    fig, ax = plt.subplots(figsize=(max(6, prop.shape[0] * 0.5 + 2), 4))
    prop.plot(kind='bar', stacked=True, ax=ax, width=0.85, legend=False)
    ax.set_ylabel('fraction of cells'); ax.set_xlabel('patient')
    ax.set_title('Per-patient LTL331-cluster assignment')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=5,
              ncol=1, frameon=False)
    savefig(fig, outdir, 'step2_per_patient_composition', plot_formats)
    plt.close(fig)
    print('[step2] wrote per-patient proportions + presence + plot')


# --------------------------------------------------------------------------- #
# Step 3 — patient-aware differential abundance
# --------------------------------------------------------------------------- #
def step3_diff_abundance(adata, args, outdir):
    assign_path = outdir / 'step1_clinical_assignments.tsv'
    if not assign_path.exists():
        print('[step3] SKIP: run step 1 first')
        return
    if not args.da_group_col or args.da_group_col not in adata.obs.columns:
        print(f'[step3] SKIP: --da-group-col not set / not in obs '
              f'(need a per-patient grouping, e.g. PRAD vs NEPC)')
        return
    df = pd.read_csv(assign_path, sep='\t', index_col=0)
    df['pred_cluster'] = df['pred_cluster'].astype(str)  # see step 2 note
    df = df[df['confidence'] >= args.conf_threshold]

    # map patient -> group
    grp = (adata.obs[[args.patient_col, args.da_group_col]]
           .drop_duplicates().set_index(args.patient_col)[args.da_group_col])
    tgt = str(args.target_cluster)
    ct = pd.crosstab(df[args.patient_col], df['pred_cluster'])
    n_tgt = ct.get(tgt, pd.Series(0, index=ct.index))
    n_tot = ct.sum(1)
    tab = pd.DataFrame({'n_target': n_tgt, 'n_total': n_tot})
    tab['group'] = grp.reindex(tab.index).values
    tab['prop'] = tab['n_target'] / tab['n_total']
    tab.to_csv(outdir / 'step3_da_input.tsv', sep='\t')

    if args.use_sccoda and have('pertpy'):
        print('[step3] pertpy available -> run scCODA externally on '
              'step3_da_input cell-counts for full compositional model; '
              'scaffold here reports the GLM + Wilcoxon fallback.')

    # Fallback A: binomial GLM with (success, failure) endog, patient = observation
    try:
        import statsmodels.api as sm
        d = tab.dropna(subset=['group'])
        d = d[d['n_total'] > 0].copy()
        if d['group'].nunique() == 2:
            d['grp_bin'] = (d['group'] == sorted(d['group'].unique())[1]).astype(int)
            endog = np.c_[d['n_target'].values,
                          (d['n_total'] - d['n_target']).values]
            X = sm.add_constant(d['grp_bin'].values.astype(float))
            res = sm.GLM(endog, X, family=sm.families.Binomial()).fit()
            with open(outdir / 'step3_glm_summary.txt', 'w') as fh:
                fh.write(res.summary().as_text())
            print(f'[step3] binomial GLM group coef p={res.pvalues[1]:.3g} '
                  f'(see step3_glm_summary.txt)')
    except Exception as e:
        print(f'[step3] GLM fallback skipped: {e}')

    # Fallback B: Mann-Whitney on per-patient proportions
    d = tab.dropna(subset=['group'])
    groups = d['group'].unique()
    if len(groups) == 2:
        a = d.loc[d['group'] == groups[0], 'prop']
        b = d.loc[d['group'] == groups[1], 'prop']
        if len(a) and len(b):
            u, p = stats.mannwhitneyu(a, b, alternative='two-sided')
            print(f'[step3] Mann-Whitney {groups[0]} vs {groups[1]} '
                  f'cluster-{tgt} proportion: U={u:.1f} p={p:.3g}')
    print('[step3] wrote DA input + tests (patient-level)')


# --------------------------------------------------------------------------- #
# Step 4 — pseudo-bulk signature scoring + null
# --------------------------------------------------------------------------- #
def step4_pseudobulk_signature(adata, args, outdir):
    if not args.signatures:
        print('[step4] SKIP: provide --signatures (cluster -> gene list)')
        return
    sigs = _load_signatures(args.signatures)
    if args.patient_col not in adata.obs.columns:
        print(f'[step4] SKIP: patient column {args.patient_col} missing')
        return

    # pseudo-bulk: mean expression per patient (clinical cells only)
    obs = adata.obs
    qry_mask = (obs[args.group_col].astype(str) != args.ref_label).values
    var_names = np.asarray(adata.var_names)
    gidx = {g: i for i, g in enumerate(var_names)}
    X = adata.layers[args.layer] if args.layer else adata.X

    patients = obs.loc[qry_mask, args.patient_col].astype(str)
    uniq = sorted(patients.unique())
    # build patient x gene pseudo-bulk (mean of normalized expression)
    import scipy.sparse as sp
    Xq = X[qry_mask]
    pb = np.zeros((len(uniq), len(var_names)), dtype=float)
    for i, pt in enumerate(uniq):
        m = (patients.values == pt)
        sub = Xq[m]
        pb[i] = np.asarray(sub.mean(0)).ravel() if sp.issparse(sub) else sub.mean(0)
    pb = pd.DataFrame(pb, index=uniq, columns=var_names)

    if args.use_decoupler and have('decoupler'):
        print('[step4] decoupler available -> recommended for AUCell/ssGSEA on '
              'the pseudo-bulk; scaffold uses the self-contained z-score + null.')

    rng = np.random.default_rng(0)
    rows = []
    for clust, genes in sigs.items():
        genes = [g for g in genes if g in gidx]
        if len(genes) < 3:
            continue
        # signature score per patient = mean z-scored expression of sig genes
        Z = (pb - pb.mean(0)) / (pb.std(0).replace(0, np.nan))
        sig_score = Z[genes].mean(1)
        # null: random matched-size gene sets
        null = np.zeros((args.n_perm, len(uniq)))
        allg = list(gidx)
        for b in range(args.n_perm):
            rs = rng.choice(allg, size=len(genes), replace=False)
            null[b] = Z[list(rs)].mean(1).values
        # per-patient empirical p (one-sided: sig > null)
        emp_p = (null >= sig_score.values).mean(0)
        for pt, sc, pv in zip(uniq, sig_score.values, emp_p):
            rows.append({'cluster_signature': clust, 'patient': pt,
                         'sig_score': sc, 'null_emp_p': pv})
    res = pd.DataFrame(rows)
    res.to_csv(outdir / 'step4_pseudobulk_signature.tsv', sep='\t', index=False)
    print(f'[step4] wrote pseudo-bulk signature scores + null '
          f'({res["cluster_signature"].nunique()} signatures x {len(uniq)} patients)')


def _load_signatures(path):
    p = Path(path)
    if p.suffix == '.json':
        return {k: list(v) for k, v in json.loads(p.read_text()).items()}
    # TSV: cluster<TAB>comma_separated_genes
    out = {}
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split('\t')
        out[parts[0]] = [g.strip() for g in parts[1].replace(',', ' ').split()]
    return out


# --------------------------------------------------------------------------- #
# Step 5 — per-patient MNN fraction (co-localization)
# --------------------------------------------------------------------------- #
def step5_mnn(adata, args, outdir):
    from sklearn.neighbors import NearestNeighbors
    if args.latent_key not in adata.obsm:
        print(f'[step5] SKIP: latent key {args.latent_key} missing')
        return
    obs = adata.obs
    Z = np.asarray(adata.obsm[args.latent_key])
    tgt = str(args.target_cluster)

    ref_mask = (obs[args.group_col].astype(str) == args.ref_label).values
    tgt_mask = ref_mask & (obs[args.ref_cluster_col].astype(str) == tgt).values
    if tgt_mask.sum() == 0:
        print(f'[step5] SKIP: no reference cluster-{tgt} cells')
        return
    Z_tgt = Z[tgt_mask]

    results = []
    clinical = obs[args.group_col].astype(str) != args.ref_label
    if args.patient_col in obs.columns:
        groups = obs.loc[clinical, args.patient_col].astype(str)
        units = [(pt, (clinical.values) & (obs[args.patient_col].astype(str) == pt).values)
                 for pt in sorted(groups.unique())]
    else:
        units = [(g, (obs[args.group_col].astype(str) == g).values)
                 for g in obs[args.group_col].unique() if g != args.ref_label]

    # kNN of target cells into the full space (for mutual check)
    nn_all = NearestNeighbors(n_neighbors=args.knn_k).fit(Z)
    for name, cmask in units:
        if cmask.sum() == 0:
            continue
        Z_c = Z[cmask]
        # target -> clinical-unit neighbors
        nn_c = NearestNeighbors(n_neighbors=args.knn_k).fit(Z_c)
        _, t2c = nn_c.kneighbors(Z_tgt)         # each target's nearest unit cells
        # clinical-unit -> global neighbors, check if a target cell is returned
        nn_t = NearestNeighbors(n_neighbors=args.knn_k).fit(Z_tgt)
        _, c2t = nn_t.kneighbors(Z_c)           # each unit cell's nearest targets
        # mutual: target i and unit j are MNN if j in t2c[i] and i in c2t[j]
        mnn_targets = 0
        c2t_sets = [set(row) for row in c2t]
        for i, neigh in enumerate(t2c):
            if any(i in c2t_sets[j] for j in neigh):
                mnn_targets += 1
        frac = mnn_targets / Z_tgt.shape[0]
        results.append({'unit': name, 'n_unit_cells': int(cmask.sum()),
                        'frac_target_with_MNN': round(frac, 4)})
    pd.DataFrame(results).to_csv(outdir / 'step5_mnn_fraction.tsv',
                                 sep='\t', index=False)
    print('[step5] wrote per-unit MNN fraction for cluster ' + tgt)

    if args.use_scib and have('scib'):
        print('[step5] scib available -> compute iLISI restricted to the '
              'cluster-17 neighborhood for the mixing-where-it-should panel.')


# --------------------------------------------------------------------------- #
# Step 6 — disease-axis stratified analysis (optional; needs --patient-metadata)
# --------------------------------------------------------------------------- #
def step6_disease_axis(args, outdir, plot_formats):
    """Stratify per-patient signature scores along a clinical disease axis.

    Optional analysis activated by --patient-metadata. Reads:
      step4_pseudobulk_signature.tsv    per-patient signature scores
      step1_clinical_assignments.tsv    per-cell label transfer (if present)
      args.patient_metadata             TSV with patient/disease_state/disease_axis
    Produces:
      step6_disease_axis_signatures.tsv     patient x signature, annotated with disease
      step6_trend_test.tsv                  Spearman + Kendall trend along disease_axis
      step6_per_patient_cluster_counts.tsv  cluster assignments per patient (annotated)
      step6_disease_heatmap.{pdf,png}       patient x signature heatmap, rows ordered
                                            by disease_axis, disease-state separators

    Gracefully no-ops with an informative message if --patient-metadata is
    missing or the required upstream outputs aren't present.
    """
    if not args.patient_metadata:
        print('[step6] SKIP: --patient-metadata not provided (this step is '
              'optional; rerun with --patient-metadata once the table is ready)')
        return

    meta_path = Path(args.patient_metadata)
    if not meta_path.exists():
        print(f'[step6] SKIP: metadata file not found: {meta_path}')
        return

    sig_path = outdir / 'step4_pseudobulk_signature.tsv'
    if not sig_path.exists():
        print(f'[step6] SKIP: {sig_path} not found (run step 4 first)')
        return

    sep = '\t' if str(meta_path).endswith(('.tsv', '.txt')) else ','
    meta = pd.read_csv(meta_path, sep=sep)
    print(f'[step6] loaded metadata: {len(meta)} rows; '
          f'columns: {list(meta.columns)}')

    pat_col = args.patient_col if args.patient_col in meta.columns else meta.columns[0]
    meta = meta.set_index(pat_col)

    missing = [c for c in (args.disease_col, args.disease_axis_col)
               if c not in meta.columns]
    if missing:
        print(f'[step6] SKIP: metadata missing required columns: {missing} '
              f'(have: {list(meta.columns)})')
        return

    meta[args.disease_axis_col] = pd.to_numeric(
        meta[args.disease_axis_col], errors='coerce')

    sig = pd.read_csv(sig_path, sep='\t')
    needed = {'patient', 'cluster_signature', 'sig_score'}
    if not needed.issubset(sig.columns):
        print(f'[step6] SKIP: step4 output missing columns {needed - set(sig.columns)}')
        return
    sig_mat = sig.pivot(index='patient', columns='cluster_signature',
                        values='sig_score')

    annot = meta[[args.disease_col, args.disease_axis_col]]
    sig_with = sig_mat.join(annot, how='left')
    n_unjoined = sig_with[args.disease_axis_col].isna().sum()
    if n_unjoined:
        print(f'[step6] {n_unjoined} patient(s) in step4 not found in '
              f'metadata; dropping them from the trend analysis')
    sig_with = sig_with.dropna(subset=[args.disease_axis_col])
    sig_with = sig_with.sort_values(args.disease_axis_col)

    if sig_with.empty:
        print('[step6] SKIP: no patients overlap between metadata and step4 output')
        return

    sig_with.to_csv(outdir / 'step6_disease_axis_signatures.tsv', sep='\t')
    print(f'[step6] joined {len(sig_with)} patients x '
          f'{sig_with.shape[1] - 2} signatures')

    # --- per-signature trend along the disease axis (Spearman + Kendall) ---
    sig_cols = [c for c in sig_with.columns
                if c not in (args.disease_col, args.disease_axis_col)]
    axis = sig_with[args.disease_axis_col].values.astype(float)
    rows = []
    for s in sig_cols:
        y = sig_with[s].values.astype(float)
        valid = ~np.isnan(y) & ~np.isnan(axis)
        if valid.sum() < 3:
            continue
        sr, sp = stats.spearmanr(axis[valid], y[valid])
        kt, kp = stats.kendalltau(axis[valid], y[valid])
        rows.append({
            'signature': s, 'n_patients': int(valid.sum()),
            'spearman_r': float(sr), 'spearman_p': float(sp),
            'kendall_tau': float(kt), 'kendall_p': float(kp),
        })
    trend = pd.DataFrame(rows).sort_values('spearman_r', ascending=False)
    trend.to_csv(outdir / 'step6_trend_test.tsv', sep='\t', index=False)
    print('[step6] disease-axis trend per signature (sorted by Spearman r):')
    print(trend.round(3).to_string(index=False))
    print('  (low n ~10 patients -> emphasise effect size & direction, '
          'not the p-values)')

    # --- per-patient cluster counts annotated with disease state ---
    asn_path = outdir / 'step1_clinical_assignments.tsv'
    if asn_path.exists():
        asn = pd.read_csv(asn_path, sep='\t', index_col=0)
        if args.patient_col in asn.columns and 'pred_cluster' in asn.columns:
            asn = asn[asn['confidence'] >= args.conf_threshold]
            counts = pd.crosstab(asn[args.patient_col],
                                 asn['pred_cluster'].astype(str))
            counts = counts.join(annot, how='left')
            if args.disease_axis_col in counts.columns:
                counts = counts.sort_values(args.disease_axis_col)
            counts.to_csv(outdir / 'step6_per_patient_cluster_counts.tsv', sep='\t')
            print(f'[step6] wrote per-patient cluster counts '
                  f'(label transfer conf >= {args.conf_threshold})')

    # --- heatmap: patients (ordered by disease_axis) x signatures (z-scored) ---
    plt = setup_style()
    M = sig_with[sig_cols].astype(float)
    mu, sd = M.mean(axis=0), M.std(axis=0).replace(0, np.nan)
    Z = (M.sub(mu, axis=1).div(sd, axis=1)).fillna(0)

    n_rows, n_cols = Z.shape
    fig, ax = plt.subplots(
        figsize=(max(6, n_cols * 0.5 + 2),
                 max(3.5, n_rows * 0.35 + 2.5))
    )
    im = ax.imshow(Z.values, aspect='auto', cmap='RdBu_r',
                   vmin=-2, vmax=2, interpolation='nearest')
    row_labels = [f'{p}  [{sig_with.loc[p, args.disease_col]}]'
                  for p in Z.index]
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(Z.columns, rotation=45, ha='right', fontsize=6)
    ax.set_title(f'Patient x LTL331-signature (z across patients)\n'
                 f'rows ordered by {args.disease_axis_col}')
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label('z-score', fontsize=6)

    # disease-state separator lines
    prev = None
    for i, p in enumerate(Z.index):
        ds = sig_with.loc[p, args.disease_col]
        if prev is not None and ds != prev:
            ax.axhline(i - 0.5, color='black', lw=0.8)
        prev = ds

    savefig(fig, outdir, 'step6_disease_heatmap', plot_formats)
    plt.close(fig)
    print('[step6] wrote disease-axis heatmap')


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument('--h5ad', required=True)
    p.add_argument('--outdir', default='patient_stats')
    p.add_argument('--latent-key', default='X_harmony',
                   help='obsm key of the integrated latent (e.g. X_harmony_dataset)')
    p.add_argument('--group-col', default='dataset')
    p.add_argument('--ref-label', default='LTL331',
                   help='value in --group-col identifying the reference (LTL331)')
    p.add_argument('--patient-col', default='patient')
    p.add_argument('--ref-cluster-col', default='original_cluster',
                   help='column holding LTL331 reference cluster labels')
    p.add_argument('--target-cluster', default='17')
    p.add_argument('--da-group-col', default=None,
                   help='per-patient grouping for Step 3 (e.g. PRAD vs NEPC)')
    p.add_argument('--signatures', default=None,
                   help='TSV (cluster<TAB>genes) or JSON for Step 4')
    p.add_argument('--layer', default=None,
                   help='layer with log-normalized expression for Step 4 '
                        '(default .X)')
    p.add_argument('--knn-k', type=int, default=30)
    p.add_argument('--conf-threshold', type=float, default=0.5)
    p.add_argument('--n-perm', type=int, default=1000)
    p.add_argument('--use-celltypist', action='store_true')
    p.add_argument('--use-decoupler', action='store_true')
    p.add_argument('--use-sccoda', action='store_true')
    p.add_argument('--use-scib', action='store_true')
    p.add_argument('--steps', nargs='+',
                   default=['1', '2', '3', '4', '5', '6'])
    p.add_argument('--patient-metadata', default=None,
                   help='OPTIONAL TSV: patient<TAB>disease_state<TAB>disease_axis'
                        '<TAB>... ; activates Step 6 disease-axis analysis. '
                        'The full pipeline runs fine without it.')
    p.add_argument('--disease-col', default='disease_state',
                   help='column in --patient-metadata with the categorical '
                        'disease label (e.g. PIN / HSPC / PRAD / CRPC / NEPC)')
    p.add_argument('--disease-axis-col', default='disease_axis',
                   help='column in --patient-metadata with an ordinal disease '
                        'rank (numeric; used for the Spearman/Kendall trend '
                        'test and the heatmap row ordering)')
    p.add_argument('--plot-format', nargs='+', default=['pdf'],
                   choices=['pdf', 'png', 'svg'])
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_adata(args.h5ad)
    print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes')
    print(f'  group counts: {dict(adata.obs[args.group_col].value_counts())}')

    steps = set(args.steps)
    if '1' in steps:
        step1_label_transfer(adata, args, outdir)
    if '2' in steps:
        step2_proportions(adata, args, outdir, args.plot_format)
    if '3' in steps:
        step3_diff_abundance(adata, args, outdir)
    if '4' in steps:
        step4_pseudobulk_signature(adata, args, outdir)
    if '5' in steps:
        step5_mnn(adata, args, outdir)
    if '6' in steps:
        print('\n=== 6: disease-axis stratified analysis (optional metadata) ===')
        step6_disease_axis(args, outdir, args.plot_format)
    print('\nDone.')


if __name__ == '__main__':
    main()
