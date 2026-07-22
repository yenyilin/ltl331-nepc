#!/bin/bash
# 02_reference_regions.sh — establish the NE-CRE vs Ad-CRE reference regions (hg38)
# that anchor both the motif test and the signal profiles.
#
# Preferred: use Baca's curated Ne-CRE (n≈14,985) / Ad-CRE (n≈4,338) BEDs if you
# downloaded them (set BACA_NECRE/BACA_ADCRE in config.sh) — they get lifted to hg38.
# Fallback: derive them from the manifest by contrasting NEPC (fg) vs adeno (bg)
# H3K27ac/FOXA peaks with bedtools (NE-specific = fg minus bg; Ad-specific = bg minus fg).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config.sh"
MAN="${MANIFEST:-$WORK/manifest.hg38.tsv}"
REF="$WORK/refregions"; mkdir -p "$REF"
CHAIN="$CHAIN_HG19_HG38"

lift_if_hg19 () { # $1 in.bed (hg19) -> $2 out.bed (hg38)
  awk 'BEGIN{OFS="\t"} !/^#/&&NF>=3{print $1,$2,$3,($4==""?"r_"NR:$4)}' "$1" > "$1.b4"
  liftOver "$1.b4" "$CHAIN" "$2" "$2.unmapped" >/dev/null 2>&1; rm -f "$1.b4"
}

if [[ -n "${BACA_NECRE:-}" && -s "$BACA_NECRE" && -n "${BACA_ADCRE:-}" && -s "$BACA_ADCRE" ]]; then
  echo "[refregions] using Baca curated CRE sets (lifting hg19->hg38)"
  lift_if_hg19 "$BACA_NECRE" "$REF/NE-CRE.hg38.bed"
  lift_if_hg19 "$BACA_ADCRE" "$REF/Ad-CRE.hg38.bed"
else
  echo "[refregions] curated CREs not provided — deriving NE-specific vs Ad-specific"
  echo "             from NEPC(fg) vs adeno(bg) H3K27ac/FOXA peaks in the manifest"
  # pool NEPC foreground peaks and adeno background peaks (prefer H3K27ac, else FOXA*)
  fg="$REF/_fg.bed"; bg="$REF/_bg.bed"; : > "$fg"; : > "$bg"
  awk -F'\t' 'NR>1 && $3=="peak" && $6=="fg"  && ($4=="H3K27ac"||$4 ~ /FOXA/){print $7}' "$MAN" | \
     while read -r f; do cut -f1-3 "$f" >> "$fg"; done
  awk -F'\t' 'NR>1 && $3=="peak" && $6=="bg"  && ($4=="H3K27ac"||$4 ~ /FOXA/){print $7}' "$MAN" | \
     while read -r f; do cut -f1-3 "$f" >> "$bg"; done
  for x in "$fg" "$bg"; do sort -k1,1 -k2,2n "$x" | bedtools merge -i - > "$x.m"; mv "$x.m" "$x"; done
  # NE-specific = fg not overlapping bg ; Ad-specific = bg not overlapping fg
  bedtools intersect -v -a "$fg" -b "$bg" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"NE_"NR}' > "$REF/NE-CRE.hg38.bed"
  bedtools intersect -v -a "$bg" -b "$fg" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Ad_"NR}' > "$REF/Ad-CRE.hg38.bed"
  rm -f "$fg" "$bg"
fi

echo "[ok] NE-CRE: $(wc -l < "$REF/NE-CRE.hg38.bed")  Ad-CRE: $(wc -l < "$REF/Ad-CRE.hg38.bed")"
echo "     -> $REF/{NE-CRE,Ad-CRE}.hg38.bed"
