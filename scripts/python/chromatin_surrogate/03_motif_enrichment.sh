#!/bin/bash
# 03_motif_enrichment.sh — THE core surrogate test: are our regulon TF motifs
# (ASCL1/NKX2-1/NEUROD1/FOXA2/...) enriched in NEPC regulatory regions vs adeno?
#
# Runs HOMER findMotifsGenome.pl (same tool Baca + Formaggio used) on the NE-CRE
# regions with Ad-CRE as background. Reports BOTH known-motif enrichment (parsed
# in step 05 for our TF set) and de-novo motifs. Optionally also runs per-dataset
# NEPC-vs-adeno FOXA peak contrasts so the result isn't anchored on one region set.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config.sh"
MAN="${MANIFEST:-$WORK/manifest.hg38.tsv}"
REF="$WORK/refregions"
MOUT="$OUT/motif"; mkdir -p "$MOUT"

run_homer () { # $1 fg.bed  $2 bg.bed  $3 label
  local fg="$1" bg="$2" label="$3" od="$MOUT/$3"
  if [[ -d "$od" && -s "$od/knownResults.txt" ]]; then echo "[skip] $label (done)"; return; fi
  echo "[homer] $label : fg=$(wc -l < "$fg") bg=$(wc -l < "$bg")"
  findMotifsGenome.pl "$fg" "$HOMER_GENOME" "$od" \
      -bg "$bg" -size 200 -mask -p "$THREADS" -nomotif >/dev/null 2>&1 || {
        # -nomotif = known-only (fast). Drop it (below) to also get de-novo motifs.
        echo "[warn] $label known-motif run failed; check HOMER genome install"; return; }
}

# 1) primary: NE-CRE (fg) vs Ad-CRE (bg)
run_homer "$REF/NE-CRE.hg38.bed" "$REF/Ad-CRE.hg38.bed" "NECRE_vs_AdCRE"

# 2) per-dataset FOXA2 NEPC-vs-adeno contrasts (robustness; only if both sides exist)
for ds in baca wang formaggio; do
  for fac in FOXA2 FOXA1; do
    fg="$WORK/_fg_${ds}_${fac}.bed"; bg="$WORK/_bg_${ds}_${fac}.bed"
    awk -F'\t' -v d="$ds" -v f="$fac" 'NR>1&&$1==d&&$3=="peak"&&$4==f&&$6=="fg"{print $7}' "$MAN" | \
       while read -r p; do cut -f1-3 "$p"; done | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null > "$fg" || true
    awk -F'\t' -v d="$ds" -v f="$fac" 'NR>1&&$1==d&&$3=="peak"&&$4==f&&$6=="bg"{print $7}' "$MAN" | \
       while read -r p; do cut -f1-3 "$p"; done | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null > "$bg" || true
    if [[ -s "$fg" && -s "$bg" ]]; then
      awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"pk_"NR}' "$fg" > "${fg}.4"; mv "${fg}.4" "$fg"
      awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"pk_"NR}' "$bg" > "${bg}.4"; mv "${bg}.4" "$bg"
      run_homer "$fg" "$bg" "${ds}_${fac}_NEPC_vs_adeno"
    fi
    rm -f "$fg" "$bg"
  done
done

# 3) cluster-17-bridge sub-analysis: Wang patient-matched DNPC vs ARPC H3K27ac.
#    201.2 (DNPC, lung met) vs 201.1 (ARPC, dura met) from one patient (Wang/Cai 2024,
#    GSE232555). DNPC = AR-/NE- mesenchymal/basal state = the human chromatin analog of the
#    cluster-17 bridge. Expect FOXA2 / AP-1 (JUN/FOSL1) / EMT (SNAI2/ZEB) motifs to enrich
#    — i.e. the cl17 program — NOT the terminal-NE drivers (kept to contrast #1). Same-patient
#    background controls for germline; the two are divergent mets, so this is a cross-sectional
#    program comparison, not an induced transition.
dfg="$WORK/_dnpc_fg.bed"; dbg="$WORK/_arpc_bg.bed"
awk -F'\t' 'NR>1 && $3=="peak" && $6=="dnpc_fg"{print $7}' "$MAN" | \
   while read -r p; do cut -f1-3 "$p"; done | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null > "$dfg" || true
awk -F'\t' 'NR>1 && $3=="peak" && $6=="arpc_bg"{print $7}' "$MAN" | \
   while read -r p; do cut -f1-3 "$p"; done | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null > "$dbg" || true
if [[ -s "$dfg" && -s "$dbg" ]]; then
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"pk_"NR}' "$dfg" > "${dfg}.4"; mv "${dfg}.4" "$dfg"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"pk_"NR}' "$dbg" > "${dbg}.4"; mv "${dbg}.4" "$dbg"
  run_homer "$dfg" "$dbg" "wang_DNPC201.2_vs_ARPC201.1_H3K27ac_cl17"
else
  echo "[skip] cl17/DNPC contrast — need wang 201.2 (dnpc_fg) + 201.1 (arpc_bg) H3K27ac peaks"
fi
rm -f "$dfg" "$dbg"

echo "[ok] HOMER knownResults under $MOUT/*/knownResults.txt"
echo "     (to also get de-novo motifs, re-run without the -nomotif flag — slower)"
