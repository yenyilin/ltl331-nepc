#!/bin/bash
# 01_liftover.sh — lift hg19 datasets (baca, wang) to hg38 so everything shares a
# build with the LTL331 data and Formaggio. BED peaks via UCSC liftOver; bigWigs
# via CrossMap (handles the bedGraph round-trip + re-binning). Formaggio (hg38)
# files are passed through unchanged.
#
# Produces lifted copies under $WORK/lifted/<dataset>/ and writes manifest.hg38.tsv
# with the `path` column rewritten to the hg38 file (build set to hg38).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/config.sh"

CHAIN="$CHAIN_HG19_HG38"
if [[ ! -s "$CHAIN" ]]; then
  echo "[setup] fetching hg19->hg38 chain"
  wget -q -O "$CHAIN" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi

OUTMAN="$WORK/manifest.hg38.tsv"
echo -e "dataset\tbuild\ttype\tfactor\tcondition\tuse\tpath" > "$OUTMAN"

# header-aware read of the manifest
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r dataset build ftype factor cond use path; do
  [[ "$use" == "ignore" ]] && continue
  dest="$WORK/lifted/$dataset"; mkdir -p "$dest"
  base="$(basename "$path")"; base="${base%.gz}"
  if [[ "$build" == "hg38" ]]; then
    # pass-through (Formaggio); gunzip if needed into dest for a uniform tree
    if [[ "$path" == *.gz ]]; then gzip -dc "$path" > "$dest/$base"; out="$dest/$base"
    else out="$path"; fi
    echo -e "$dataset\thg38\t$ftype\t$factor\t$cond\t$use\t$out" >> "$OUTMAN"
    continue
  fi
  # ---- hg19 -> hg38 ----
  if [[ "$ftype" == "peak" ]]; then
    src="$dest/${base}.in"; [[ "$path" == *.gz ]] && gzip -dc "$path" > "$src" || cp "$path" "$src"
    # keep only chrom,start,end,name (col4 may be absent in broadPeak summit cols)
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1,$2,$3,($4==""?"pk_"NR:$4)}' "$src" > "${src}.bed4"
    out="$dest/${base%.*}.hg38.bed"
    liftOver "${src}.bed4" "$CHAIN" "$out" "${out}.unmapped" >/dev/null 2>&1 || true
    n_in=$(wc -l < "${src}.bed4"); n_out=$(wc -l < "$out" 2>/dev/null || echo 0)
    echo "[liftOver peak] $dataset/$base  $n_in -> $n_out"
    rm -f "$src" "${src}.bed4"
    echo -e "$dataset\thg38\t$ftype\t$factor\t$cond\t$use\t$out" >> "$OUTMAN"
  else  # signal bigWig
    src="$path"; [[ "$path" == *.gz ]] && { gzip -dc "$path" > "$dest/$base"; src="$dest/$base"; }
    prefix="$dest/${base%.*}.hg38"
    # CrossMap v0.6+: `CrossMap bigwig`; older: `CrossMap.py bigwig`
    if command -v CrossMap >/dev/null 2>&1; then CM="CrossMap"; else CM="CrossMap.py"; fi
    $CM bigwig "$CHAIN" "$src" "$prefix" >/dev/null 2>&1
    out="${prefix}.bw"
    echo "[liftOver bw]   $dataset/$base -> $(basename "$out")"
    echo -e "$dataset\thg38\t$ftype\t$factor\t$cond\t$use\t$out" >> "$OUTMAN"
  fi
done

echo "[ok] wrote $OUTMAN — all rows now hg38. Use this as MANIFEST downstream:"
echo "     export MANIFEST=$OUTMAN"
