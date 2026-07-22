#!/bin/bash
# 04_signal_profiles.sh — deepTools signal evidence: do FOXA1/FOXA2/H3K27ac/ASCL1
# ChIP tracks pile up over NE-CRE vs Ad-CRE regions? Uses the hg38 bigWigs from the
# manifest (use==signal) and the reference regions from step 02. Produces a matrix,
# a heatmap, and a profile per region set — the second chromatin-surrogate panel.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config.sh"
MAN="${MANIFEST:-$WORK/manifest.hg38.tsv}"
REF="$WORK/refregions"
SOUT="$OUT/signal"; mkdir -p "$SOUT"

# collect signal bigWigs + readable labels (dataset:factor:condition)
BWS=(); LABELS=()
while IFS=$'\t' read -r dataset build ftype factor cond use path; do
  [[ "$ftype" == "signal" && "$use" == "signal" && -s "$path" ]] || continue
  BWS+=("$path"); LABELS+=("${dataset}:${factor}:${cond}")
done < <(tail -n +2 "$MAN")

if [[ ${#BWS[@]} -eq 0 ]]; then echo "[err] no signal bigWigs in manifest"; exit 1; fi
echo "[signal] ${#BWS[@]} bigWigs over NE-CRE / Ad-CRE"

computeMatrix reference-point --referencePoint center \
  -S "${BWS[@]}" \
  -R "$REF/NE-CRE.hg38.bed" "$REF/Ad-CRE.hg38.bed" \
  -a 2000 -b 2000 --binSize 50 --missingDataAsZero --skipZeros \
  --samplesLabel "${LABELS[@]}" \
  -p "$THREADS" -o "$SOUT/matrix.gz" --outFileNameMatrix "$SOUT/matrix.tab"

plotHeatmap -m "$SOUT/matrix.gz" -o "$SOUT/heatmap_NEvsAd.pdf" \
  --regionsLabel "NE-CRE" "Ad-CRE" --colorMap RdBu_r \
  --heatmapHeight 12 --refPointLabel center --dpi 300

plotProfile -m "$SOUT/matrix.gz" -o "$SOUT/profile_NEvsAd.pdf" \
  --regionsLabel "NE-CRE" "Ad-CRE" --perGroup --numPlotsPerRow 4 --dpi 300

echo "[ok] $SOUT/{heatmap_NEvsAd.pdf,profile_NEvsAd.pdf,matrix.tab}"
