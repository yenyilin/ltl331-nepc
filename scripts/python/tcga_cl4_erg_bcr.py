#!/usr/bin/env python3
"""
tcga_cl4_erg_bcr.py - does the LTL331 cluster-4 signature stratify BCR within ERG-fusion strata?

External cohort: TCGA-PRAD (PanCancer Atlas, cBioPortal bundle).
- Expression : data_mrna_seq_v2_rsem.txt        -> cluster-4 mean-z signature score
- ERG status : data_sv.txt (TMPRSS2:ERG fusion) -> ERG+ / ERG- strata
- Endpoint   : PFS_STATUS / PFS_MONTHS          -> BCR proxy (PFI; DFS has only 30 events)

Within each ERG stratum we test cluster-4-high vs -low (KM + log-rank) and the continuous
score (Cox HR), plus a cl4 x ERG interaction term. Validation, not discovery: report HRs
with CIs and do not over-read a null.

Outputs:
  figures/figS_tcga_cl4_erg_bcr.{pdf,png}
  data/tcga_cl4_erg_bcr_stats.tsv

Paths:
  ROOT     repository root; auto-detected from this script's location (scripts/ -> root).
           Override with `export LTL331_ROOT=/path/to/repo` if run from another directory.
  SIG_FILE in-repo: data/chea3_input/cl4_top150.txt  (cluster-4 top-150 signature).
  FIGDIR   output: figures/  (created if missing).
  TCGA     TCGA-PRAD PanCancer Atlas (2018) cBioPortal bundle, NOT included in the repo
           (size / redistribution). Download it and unpack to
           public-data/tcga/prad_tcga_pan_can_atlas_2018/ :
           https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# Repo root: auto-detected from this script's location (scripts/ -> repo root).
# Override with `export LTL331_ROOT=/path/to/repo` if you run from elsewhere.
ROOT = os.environ.get(
    "LTL331_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
TCGA = os.path.join(ROOT, "public-data/tcga/prad_tcga_pan_can_atlas_2018")
SIG_FILE = os.path.join(ROOT, "data/chea3_input/cl4_top150.txt")
FIGDIR = os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)


def patient_of(barcode):
    """TCGA-XX-XXXX-01 -> TCGA-XX-XXXX"""
    return "-".join(str(barcode).split("-")[:3])


def is_primary(barcode):
    """sample type -01 = primary solid tumor"""
    parts = str(barcode).split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


# ---------------------------------------------------------------- signature
sig = [g.strip() for g in open(SIG_FILE) if g.strip()]
print(f"[sig] cluster-4 signature: {len(sig)} genes")

# ---------------------------------------------------------------- expression
print("[expr] loading RSEM matrix (this is the slow step)...")
expr = pd.read_csv(os.path.join(TCGA, "data_mrna_seq_v2_rsem.txt"),
                   sep="\t", low_memory=False)
expr = expr.drop(columns=[c for c in ("Entrez_Gene_Id",) if c in expr.columns])
expr = expr[expr["Hugo_Symbol"].notna()]
# collapse duplicate symbols by mean, genes as index
expr = expr.groupby("Hugo_Symbol").mean(numeric_only=True)
# keep primary-tumor sample columns only
prim_cols = [c for c in expr.columns if is_primary(c)]
expr = expr[prim_cols]
print(f"[expr] {expr.shape[0]} genes x {expr.shape[1]} primary-tumor samples")

sig_in = [g for g in sig if g in expr.index]
print(f"[sig] {len(sig_in)}/{len(sig)} signature genes present in TCGA expression")

# mean-z module score: log2(RSEM+1) -> per-gene z across samples -> mean over sig genes
logx = np.log2(expr.loc[sig_in] + 1.0)
z = logx.sub(logx.mean(axis=1), axis=0).div(logx.std(axis=1).replace(0, np.nan), axis=0)
cl4_score = z.mean(axis=0)                       # per sample
cl4 = pd.DataFrame({"sample": cl4_score.index, "cl4_score": cl4_score.values})
cl4["patient"] = cl4["sample"].map(patient_of)
# collapse to patient (mean if multiple primary samples)
cl4 = cl4.groupby("patient", as_index=False)["cl4_score"].mean()

# ---------------------------------------------------------------- ERG status
sv = pd.read_csv(os.path.join(TCGA, "data_sv.txt"), sep="\t", low_memory=False)
sv["patient"] = sv["Sample_Id"].map(patient_of)
sv_patients = set(sv["patient"])
erg_pos = set(sv.loc[(sv["Site1_Hugo_Symbol"] == "ERG") |
                     (sv["Site2_Hugo_Symbol"] == "ERG"), "patient"])
erg = pd.DataFrame({"patient": sorted(sv_patients)})
erg["ERG"] = np.where(erg["patient"].isin(erg_pos), "ERG+", "ERG-")
print(f"[erg] ERG+ {len(erg_pos)} / SV-called {len(sv_patients)} patients "
      f"({100*len(erg_pos)/len(sv_patients):.0f}%)")

# ---------------------------------------------------------------- clinical / PFS
clin = pd.read_csv(os.path.join(TCGA, "data_clinical_patient.txt"),
                   sep="\t", skiprows=4, low_memory=False)
clin = clin.rename(columns={"PATIENT_ID": "patient"})
clin["PFS_event"] = clin["PFS_STATUS"].astype(str).str.startswith("1").astype(int)
clin["PFS_months"] = pd.to_numeric(clin["PFS_MONTHS"], errors="coerce")
clin = clin[["patient", "PFS_event", "PFS_months"]].dropna(subset=["PFS_months"])

# ---------------------------------------------------------------- merge
df = cl4.merge(erg, on="patient", how="inner").merge(clin, on="patient", how="inner")
df = df[df["PFS_months"] > 0]
med = df["cl4_score"].median()
df["cl4_group"] = np.where(df["cl4_score"] >= med, "cl4-high", "cl4-low")
print(f"[merge] analyzable patients: {len(df)}  "
      f"(ERG+ {sum(df.ERG=='ERG+')}, ERG- {sum(df.ERG=='ERG-')}); "
      f"PFS events {int(df.PFS_event.sum())}")

# ---------------------------------------------------------------- stats
rows = []


def cox_continuous(sub, label):
    """Cox on continuous (z-scaled) cl4 score."""
    s = sub[["PFS_months", "PFS_event", "cl4_score"]].copy()
    s["cl4_score"] = (s["cl4_score"] - s["cl4_score"].mean()) / s["cl4_score"].std()
    n, ev = len(s), int(s["PFS_event"].sum())
    try:
        cph = CoxPHFitter().fit(s, "PFS_months", "PFS_event")
        r = cph.summary.loc["cl4_score"]
        return dict(stratum=label, n=n, events=ev, test="Cox continuous (per SD)",
                    HR=float(r["exp(coef)"]),
                    HR_lo=float(r["exp(coef) lower 95%"]),
                    HR_hi=float(r["exp(coef) upper 95%"]),
                    p=float(r["p"]))
    except Exception as e:
        return dict(stratum=label, n=n, events=ev, test="Cox continuous (per SD)",
                    HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=np.nan, note=str(e)[:60])


def logrank_high_low(sub, label):
    hi = sub[sub.cl4_group == "cl4-high"]
    lo = sub[sub.cl4_group == "cl4-low"]
    n, ev = len(sub), int(sub["PFS_event"].sum())
    try:
        lr = logrank_test(hi["PFS_months"], lo["PFS_months"],
                          hi["PFS_event"], lo["PFS_event"])
        return dict(stratum=label, n=n, events=ev, test="log-rank high vs low",
                    HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=float(lr.p_value))
    except Exception as e:
        return dict(stratum=label, n=n, events=ev, test="log-rank high vs low",
                    HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=np.nan, note=str(e)[:60])


for label, sub in [("ERG+", df[df.ERG == "ERG+"]),
                   ("ERG-", df[df.ERG == "ERG-"]),
                   ("all", df)]:
    rows.append(cox_continuous(sub, label))
    rows.append(logrank_high_low(sub, label))

# interaction: cl4 x ERG in a combined Cox
inter = df[["PFS_months", "PFS_event", "cl4_score"]].copy()
inter["cl4_score"] = (inter["cl4_score"] - inter["cl4_score"].mean()) / inter["cl4_score"].std()
inter["ERG_pos"] = (df["ERG"] == "ERG+").astype(int).values
inter["cl4_x_ERG"] = inter["cl4_score"] * inter["ERG_pos"]
try:
    cph = CoxPHFitter().fit(inter, "PFS_months", "PFS_event")
    ip = float(cph.summary.loc["cl4_x_ERG", "p"])
    rows.append(dict(stratum="interaction", n=len(inter), events=int(inter.PFS_event.sum()),
                     test="Cox cl4 x ERG term", HR=np.nan, HR_lo=np.nan, HR_hi=np.nan, p=ip))
except Exception as e:
    ip = np.nan
    rows.append(dict(stratum="interaction", n=len(inter), events=int(inter.PFS_event.sum()),
                     test="Cox cl4 x ERG term", HR=np.nan, HR_lo=np.nan, HR_hi=np.nan,
                     p=np.nan, note=str(e)[:60]))

stats = pd.DataFrame(rows)
stats_path = os.path.join(ROOT, "data/tcga_cl4_erg_bcr_stats.tsv")
stats.to_csv(stats_path, sep="\t", index=False)
print("\n=== STATS ===")
print(stats.to_string(index=False))
print(f"\n[out] {stats_path}")

# ---------------------------------------------------------------- KM figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, label in zip(axes, ["ERG+", "ERG-"]):
    sub = df[df.ERG == label]
    for grp, color in [("cl4-high", "#c0392b"), ("cl4-low", "#2c6fbb")]:
        g = sub[sub.cl4_group == grp]
        if len(g) == 0:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(g["PFS_months"], g["PFS_event"], label=f"{grp} (n={len(g)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    srow = stats[(stats.stratum == label) & (stats.test.str.startswith("Cox"))].iloc[0]
    lrp = stats[(stats.stratum == label) & (stats.test.str.startswith("log-rank"))].iloc[0]["p"]
    ax.set_title(f"{label}  (n={len(sub)}, events={int(sub.PFS_event.sum())})\n"
                 f"Cox HR/SD={srow['HR']:.2f} [{srow['HR_lo']:.2f}–{srow['HR_hi']:.2f}], "
                 f"p={srow['p']:.2g}; log-rank p={lrp:.2g}", fontsize=13)
    ax.set_xlabel("Months (PFS / BCR proxy)", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=11, loc="lower left")
axes[0].set_ylabel("Progression-free probability", fontsize=13)
fig.suptitle("TCGA-PRAD: cluster-4 signature vs BCR (PFS), stratified by TMPRSS2:ERG fusion "
             f"— interaction p={ip:.2g}", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIGDIR, f"figS_tcga_cl4_erg_bcr.{ext}"), dpi=200,
                bbox_inches="tight")
print(f"[out] {os.path.join(FIGDIR, 'figS_tcga_cl4_erg_bcr.pdf')}")
