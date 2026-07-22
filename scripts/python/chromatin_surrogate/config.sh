#!/bin/bash
# config.sh — central paths + parameters for the chromatin-surrogate pipeline.
# EDIT the four RAW_* paths and ROOT to match where you untarred the GEO downloads,
# then `source config.sh` (the numbered scripts source it themselves).
set -euo pipefail

# ---- where the GEO downloads live (the *_RAW.tar you already fetched) ----------
export ROOT="${ROOT:-/path/to/chromatin_surrogate}"     # working root on the node
export RAW_BACA="${RAW_BACA:-$ROOT/raw/GSE161948}"      # untarred GSE161948_RAW.tar (hg19; BED+BW)
export RAW_WANG="${RAW_WANG:-$ROOT/raw/GSE232555}"      # untarred GSE232555_RAW.tar (hg19; BED+BW)
export RAW_FORM="${RAW_FORM:-$ROOT/raw/GSE281784}"      # untarred GSE281784_RAW.tar (hg38; BW only)
# Baca curated regulatory elements (from the paper's Supplementary Data / GEO supp).
# If you have them, point here (hg19 BED). If empty, 02_reference_regions.sh derives them.
export BACA_NECRE="${BACA_NECRE:-}"                      # e.g. $ROOT/raw/GSE161948/Ne-CRE.bed  (n≈14,985)
export BACA_ADCRE="${BACA_ADCRE:-}"                      # e.g. $ROOT/raw/GSE161948/Ad-CRE.bed  (n≈4,338)

# ---- genome -------------------------------------------------------------------
export TARGET_BUILD="hg38"                               # LTL331 + Formaggio are GRCh38
export HOMER_GENOME="hg38"                               # must be installed: see README step 0
export CHAIN_HG19_HG38="${CHAIN_HG19_HG38:-$ROOT/ref/hg19ToHg38.over.chain.gz}"
export HG38_FA="${HG38_FA:-$ROOT/ref/hg38.fa}"           # only needed for AME/FIMO (optional)

# ---- outputs ------------------------------------------------------------------
export WORK="${WORK:-$ROOT/work}"                        # organized + lifted files
export OUT="${OUT:-$ROOT/out}"                           # motif/signal/figure outputs
export MANIFEST="${MANIFEST:-$WORK/manifest.tsv}"
export THREADS="${THREADS:-16}"

# ---- our regulon TF set (the question: do these gain accessibility/binding in NEPC?) --
# Used by 05_make_s15b_chromatin.py to pull the matching rows out of HOMER results.
export REGULON_TFS="MSX1 ASCL1 NEUROD1 FOXA2 FOXA1 NKX2-1 POU3F2 INSM1 SOX2 ONECUT2 REST RUNX3 SNAI2 TWIST1"

mkdir -p "$WORK" "$OUT" "$ROOT/ref"
echo "[config] ROOT=$ROOT TARGET_BUILD=$TARGET_BUILD THREADS=$THREADS"
