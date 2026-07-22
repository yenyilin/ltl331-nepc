#!/bin/bash
# run_all.sh — one-command, fully-provenanced entry point for the chromatin-surrogate
# analysis (Supplementary Fig. S15B). This IS the workflow: linear, modular, idempotent,
# and it records everything a reviewer needs to reproduce the run — code commit, tool
# versions, input checksums, parameters, seed, environment lock. No workflow engine
# required; the provenance it writes is the reproducibility artifact.
#
# Usage:
#   bash run_all.sh                # stops after 00 for manifest review (recommended)
#   bash run_all.sh --yes          # skip the manual-review pause (manifest already vetted)
#   bash run_all.sh --from 02      # resume at a given step (idempotent; safe to re-run)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source "$HERE/config.sh"

ASSUME_YES=0; FROM=00
while [[ $# -gt 0 ]]; do case "$1" in
  --yes) ASSUME_YES=1;; --from) FROM="$2"; shift;; *) echo "unknown arg $1"; exit 2;; esac; shift; done

export SEED="${SEED:-0}"                    # pin any stochastic step
LOG="$OUT/run_$(date +%Y%m%d_%H%M%S).log"
PROV="$OUT/provenance.txt"
mkdir -p "$OUT"

# ---- provenance capture (the part that prevents "is this reproducible?" pushback) ----
{
  echo "# provenance — chromatin-surrogate pipeline"
  echo "date:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host:        $(hostname)"
  echo "user:        $(whoami)"
  echo "seed:        $SEED"
  echo "target_build:$TARGET_BUILD"
  echo "code_commit: $(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo 'not-a-git-checkout')"
  echo "code_status: $(git -C "$HERE" status --porcelain 2>/dev/null | wc -l) uncommitted file(s)"
  echo "## tool versions"
  for t in liftOver CrossMap bedtools computeMatrix findMotifsGenome.pl python3; do
    v=$("$t" --version 2>&1 | head -1 || echo 'n/a'); echo "  $t: $v"
  done
  echo "## conda env lock -> env.lock.yml"
  conda env export 2>/dev/null > "$HERE/env.lock.yml" && echo "  written env.lock.yml" || echo "  (conda export unavailable)"
  echo "## input checksums (raw GEO files)"
  for d in "$RAW_BACA" "$RAW_WANG" "$RAW_FORM"; do
    [[ -d "$d" ]] && find "$d" -type f \( -name '*.bed*' -o -name '*.bw' -o -name '*Peak*' -o -name '*.bigWig' \) \
        -exec md5sum {} \; 2>/dev/null | sort
  done
} | tee "$PROV"
echo "[provenance] -> $PROV"; echo

# ---- the pipeline (idempotent; each step skips if its outputs exist) ----
step_ge () { [[ "$FROM" -le "$1" ]]; }   # numeric compare (00,01,...)

exec > >(tee -a "$LOG") 2>&1   # everything from here is logged
echo "=== run_all START $(date) | log=$LOG ==="

if step_ge 00; then echo ">>> 00 organize"; python3 00_organize.py; fi

if [[ $ASSUME_YES -eq 0 && "$FROM" -le 00 ]]; then
  echo
  echo ">>> PAUSE: review $MANIFEST — fix condition=UNKNOWN / use columns"
  echo "    (Baca LuCaP NEPC-vs-adeno tags need a manual lookup), then re-run with:"
  echo "    bash run_all.sh --yes --from 01"
  exit 0
fi

if step_ge 01; then echo ">>> 01 liftover";          bash 01_liftover.sh; fi
export MANIFEST="$WORK/manifest.hg38.tsv"
if step_ge 02; then echo ">>> 02 reference regions"; bash 02_reference_regions.sh; fi
if step_ge 03; then echo ">>> 03 motif enrichment";  bash 03_motif_enrichment.sh; fi
if step_ge 04; then echo ">>> 04 signal profiles";   bash 04_signal_profiles.sh; fi
if step_ge 05; then echo ">>> 05 figure";            python3 05_make_s15b_chromatin.py --outdir "$OUT/figures"; fi

echo "=== run_all DONE $(date) ==="
echo "figure:     $OUT/figures/s15b_chromatin_motif.{pdf,png}"
echo "provenance: $PROV   env lock: $HERE/env.lock.yml   log: $LOG"
