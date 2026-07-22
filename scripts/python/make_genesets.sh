#!/usr/bin/env bash
#
# make_genesets.sh — build the gene-set TSVs for pathway_revision.py (R2 #5 + cluster 6).
# Pins MSigDB v7.5.1 (Methods/reproducibility). Two outputs:
#   1) genesets_hallmark50.tsv  — all 50 Hallmark sets (HALLMARK_ prefix stripped),
#                                 drives the per-cluster Hallmark heatmap (Fig 3D reproduction)
#   2) genesets_cluster6.tsv    — focused panel for the per-timepoint dynamics
#                                 (clusters 3/4/6/12): 15 relevant Hallmark + 4 cilium sets
#
# GMT format is  name<TAB>description<TAB>gene1<TAB>...  — we drop the description and
# emit the script's expected  name<TAB>comma-separated-genes.  HALLMARK_ prefix is
# stripped on output for clean labels; GO names (GOBP_/GOCC_) pass through as-is.
#
# Usage:
#   make_genesets.sh /path/to/msigdb_v7.5.1_GMTs
#
set -euo pipefail
GMT="${1:?usage: make_genesets.sh <dir with *.v7.5.1.symbols.gmt>}"
H="$GMT/h.all.v7.5.1.symbols.gmt"
# C5 GO: prefer the combined c5.all (has GOBP+GOCC); fall back to split files.
if [ -f "$GMT/c5.all.v7.5.1.symbols.gmt" ]; then
  C5=("$GMT/c5.all.v7.5.1.symbols.gmt")
else
  C5=("$GMT/c5.go.cc.v7.5.1.symbols.gmt" "$GMT/c5.go.bp.v7.5.1.symbols.gmt")
fi

# extract_named <wanted_file> <strip_prefix:yes|no> <gmt...> -> stdout (name<TAB>genes)
extract_named() {
  local want="$1" strip="$2"; shift 2
  awk -F'\t' -v strip="$strip" '
    NR==FNR { sub(/\r$/,""); if($0!~/^#/ && $0!="") keep[$0]=1; next }
    { sub(/\r$/,"") }
    keep[$1] {
      name=$1; if(strip=="yes") sub(/^HALLMARK_/,"",name)
      printf "%s\t", name
      for(i=3;i<=NF;i++) printf "%s%s",(i>3?",":""),$i
      print ""
    }' "$want" "$@"
}

# --- 1) all 50 Hallmark (strip prefix) — no wanted-list needed, take the whole file ---
awk -F'\t' '{sub(/\r$/,""); name=$1; sub(/^HALLMARK_/,"",name);
  printf "%s\t",name; for(i=3;i<=NF;i++)printf "%s%s",(i>3?",":""),$i; print ""}' \
  "$H" > genesets_hallmark50.tsv
echo "genesets_hallmark50.tsv : $(wc -l < genesets_hallmark50.tsv) sets (expect 50)"

# --- 2) focused panel for dynamics (clusters 3/4/6/12) ---
cat > /tmp/want_hallmark.$$ <<'EOF'
HALLMARK_KRAS_SIGNALING_UP
HALLMARK_KRAS_SIGNALING_DN
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
HALLMARK_ANGIOGENESIS
HALLMARK_TNFA_SIGNALING_VIA_NFKB
HALLMARK_HYPOXIA
HALLMARK_IL6_JAK_STAT3_SIGNALING
HALLMARK_MYC_TARGETS_V1
HALLMARK_MYC_TARGETS_V2
HALLMARK_E2F_TARGETS
HALLMARK_G2M_CHECKPOINT
HALLMARK_ANDROGEN_RESPONSE
HALLMARK_PANCREAS_BETA_CELLS
HALLMARK_P53_PATHWAY
HALLMARK_APOPTOSIS
EOF
cat > /tmp/want_cilium.$$ <<'EOF'
GOCC_CILIUM
GOCC_NON_MOTILE_CILIUM
GOCC_MOTILE_CILIUM
GOBP_CILIUM_ORGANIZATION
EOF

{
  extract_named /tmp/want_hallmark.$$ yes "$H"
  extract_named /tmp/want_cilium.$$   no  "${C5[@]}"
} > genesets_cluster6.tsv
rm -f /tmp/want_hallmark.$$ /tmp/want_cilium.$$
echo "genesets_cluster6.tsv   : $(wc -l < genesets_cluster6.tsv) sets (expect 19)"
echo "  names:"; cut -f1 genesets_cluster6.tsv | sed 's/^/    /'

# Sanity: any wanted set that didn't match (usually a GO-name typo)
echo ">> if a count is short, a name didn't match the GMT — re-grep the C5 file."
