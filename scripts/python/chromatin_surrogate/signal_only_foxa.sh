#!/bin/bash
# signal_only_foxa.sh — LEAN supplementary panel: FOXA2/ASCL1/FOXA1 ChIP signal over
# NE-CRE vs Ad-CRE regions. The targeted re-analysis that covers what the Shrestha
# clinical footprint set LACKS (FOXA2, SOX2, direct ASCL1 binding) — without the full
# 00-05 pipeline or its manifest curation. NO HOMER motif step (Shrestha Supp S5 covers
# motifs). deepTools signal only. ~1-2 h on the node.
#
# Inputs (env or config.sh): a track manifest (signal_foxa.tracks.tsv), hg38 NE-CRE /
# Ad-CRE BEDs, hg19->hg38 chain. hg19 tracks/regions are lifted automatically.
#
# Usage:
#   source config.sh                       # for CHAIN/WORK/OUT/THREADS (or export them)
#   export NECRE=.../NE-CRE.bed ADCRE=.../Ad-CRE.bed   # hg38 preferred; hg19 auto-lifted
#   bash signal_only_foxa.sh signal_foxa.tracks.tsv
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$HERE/config.sh" ]] && source "$HERE/config.sh" || true
TRACKS="${1:-$HERE/signal_foxa.tracks.tsv}"
OUT="${OUT:-$HERE/out}"; WORK="${WORK:-$HERE/work}"; THREADS="${THREADS:-8}"
CHAIN="${CHAIN_HG19_HG38:-$HERE/ref/hg19ToHg38.over.chain.gz}"
SOUT="$OUT/signal_foxa"; LIFT="$WORK/lifted_foxa"; mkdir -p "$SOUT" "$LIFT"
CM=$(command -v CrossMap >/dev/null 2>&1 && echo CrossMap || echo CrossMap.py)

# ---- reference regions (lift to hg38 if needed) ----
[[ -n "${NECRE:-}" && -s "$NECRE" ]] || { echo "[err] set NECRE=<NE-CRE bed>"; exit 2; }
[[ -n "${ADCRE:-}" && -s "$ADCRE" ]] || { echo "[err] set ADCRE=<Ad-CRE bed>"; exit 2; }
lift_bed () { # $1 in  $2 out  $3 build
  if [[ "$3" == "hg19" ]]; then
    awk 'BEGIN{OFS="\t"}!/^#/&&NF>=3{print $1,$2,$3,"r_"NR}' "$1" > "$2.in"
    liftOver "$2.in" "$CHAIN" "$2" "$2.unmapped" >/dev/null 2>&1; rm -f "$2.in"
  else cp "$1" "$2"; fi
}
NECRE_B="${NECRE_BUILD:-hg38}"; ADCRE_B="${ADCRE_BUILD:-hg38}"
lift_bed "$NECRE" "$LIFT/NE-CRE.hg38.bed" "$NECRE_B"
lift_bed "$ADCRE" "$LIFT/Ad-CRE.hg38.bed" "$ADCRE_B"

# ---- bigWigs: lift hg19 -> hg38, collect hg38 paths + labels ----
BWS=(); LABELS=()
while IFS=$'\t' read -r path label build; do
  [[ "$path" =~ ^#|^$ ]] && continue
  [[ -s "$path" ]] || { echo "[warn] missing track: $path (skipped)"; continue; }
  if [[ "$build" == "hg19" ]]; then
    out="$LIFT/$(basename "${path%.*}").hg38"
    echo "[lift] $label"; $CM bigwig "$CHAIN" "$path" "$out" >/dev/null 2>&1
    BWS+=("$out.bw")
  else
    BWS+=("$path")
  fi
  LABELS+=("$label")
done < "$TRACKS"
[[ ${#BWS[@]} -gt 0 ]] || { echo "[err] no usable tracks in $TRACKS"; exit 1; }
echo "[signal] ${#BWS[@]} tracks: ${LABELS[*]}"

# ---- provenance ----
{ echo "# signal_only_foxa provenance  $(date -u +%FT%TZ)"
  echo "tracks: $TRACKS"; echo "NECRE: $NECRE ($NECRE_B)  ADCRE: $ADCRE ($ADCRE_B)"
  for t in "${BWS[@]}"; do md5sum "$t" 2>/dev/null; done
} > "$SOUT/provenance.txt"

# ---- deepTools ----
computeMatrix reference-point --referencePoint center \
  -S "${BWS[@]}" -R "$LIFT/NE-CRE.hg38.bed" "$LIFT/Ad-CRE.hg38.bed" \
  -a 2000 -b 2000 --binSize 50 --missingDataAsZero --skipZeros \
  --samplesLabel "${LABELS[@]}" -p "$THREADS" -o "$SOUT/matrix.gz"

plotHeatmap -m "$SOUT/matrix.gz" -o "$SOUT/foxa_signal_heatmap.pdf" \
  --regionsLabel "NE-CRE" "Ad-CRE" --colorMap Reds --heatmapHeight 12 \
  --refPointLabel center --dpi 300
plotProfile -m "$SOUT/matrix.gz" -o "$SOUT/foxa_signal_profile.pdf" \
  --regionsLabel "NE-CRE" "Ad-CRE" --perGroup --numPlotsPerRow 3 --dpi 300

echo "[ok] $SOUT/{foxa_signal_heatmap.pdf,foxa_signal_profile.pdf}  + provenance.txt"
echo "     (supplementary: FOXA2/ASCL1/FOXA1 signal enriched at NE-CRE vs Ad-CRE)"
