# Chromatin-surrogate pipeline (Supplementary Fig. S15B)

**Goal.** No matched ATAC-seq is available for the LTL331 model. This pipeline asks, from
public bulk ChIP/ATAC in prostate NEPC vs adenocarcinoma, whether the nominated regulon TFs
(ASCL1, NKX2-1, NEUROD1, FOXA2, MSX1, POU3F2, INSM1, SOX2, ONECUT2, REST, …) show
**concordant motif enrichment and ChIP signal in NEPC regulatory regions**. This provides
chromatin-*consistent* (not chromatin-*proven*) evidence.

**Datasets (public, from GEO):**
| dataset | accession | build | files | role |
|---|---|---|---|---|
| Baca 2021 | GSE161948 | **hg19** (confirmed) | BED + bigWig: **FOXA1, H3K27ac, H3K4me3, H3K27me3** (29 LuCaP lines; **no ASCL1, no ATAC** in the RAW.tar) | **primary peak/motif anchor**; NE-CRE/Ad-CRE from NEPC-vs-adeno H3K27ac/FOXA1. NE-vs-adeno is the GEO `histology:` field (NE = 49, 93, 145.1, 145.2, 173.1, 208.1; all other LuCaP = adeno — **note 173.2 = adeno, not NE**) |
| Wang/Cai 2024 | GSE232555 | **hg19** (confirmed: chr1=249,250,621) | BED + bigWig: FOXA2, FOXA1, AR, JUN/FOSL1, H3K27ac, H3K4me2, ATAC | **mechanistic FOXA2/AP-1 study, NOT a cohort** — baseline arms only (drop all si*/OE/KO/drug). Subtypes: **H660 = NEPC**; **PC3 + PDX-201.2 = DNPC** (AR−/NE−, the human analog of cl17); **PDX-201.1 + LNCaP = ARPC**. Contributes (i) signal over Baca NE-CRE; (ii) the cl17-bridge contrast 201.2(DNPC) vs 201.1(ARPC) H3K27ac (patient-matched mets) |
| Formaggio 2025 | GSE281784 | **hg38** | **bigWig only** | FOXA1/FOXA2/H3K27ac signal in NEPC (H660, LuCaP145 PDX, MSKPCa1) |

Wang is **GSM-allow-listed** in `00_organize.py` (`WANG_GSM`), not filename-tagged — its bare-numeric
PDX names + perturbation arms make filename heuristics unsafe. Non-allow-listed GSMs → `IGNORE`. So
`wang` has **fg=bg=0 by design**; it appears only as `signal` + the `dnpc_fg`/`arpc_bg` matched pair.

`GSE281784` is the Formaggio **ChIP SubSeries** (matches the untarred dir); the paper prose cites the
SuperSeries `GSE281794` — one digit apart, don't "correct" the config. ASCL1 ChIP is **not** in any of
these tars (use GSE183200/GSE156290, both hg38, only if one wants direct ASCL1 tracks).

Build differs (hg19 vs hg38) → **step 01 lifts the hg19 sets to hg38** (= LTL331 build); Formaggio
(hg38) passes through. Baca and Wang are **both confirmed hg19** (Wang verified empirically:
`bigWigInfo` chr1 = 249,250,621; basesCovered = 3,137,161,264; chromCount 93). For any NEW bigWig of
unknown build, still check before lifting — lifting an hg38 file corrupts coordinates silently — via
`bigWigInfo <file>.bw` (chrom sizes) and fix the `build` column in `00_organize.py`'s `DATASETS` if wrong.
Formaggio has no peaks, so its bigWigs are *signal* only; peaks/motifs come from Baca + Wang.

## Environment (once)
```bash
conda env create -f env.yml && conda activate chromsurr
# HOMER genome is a separate download (~1 GB), not in conda:
perl $(dirname $(which findMotifsGenome.pl))/../share/homer*/configureHomer.pl -install hg38
```

## Configure (once)
Edit `config.sh`: set `ROOT`, and the `RAW_BACA / RAW_WANG / RAW_FORM` dirs to where
you untarred each `*_RAW.tar`. If you also grabbed Baca's curated **Ne-CRE / Ad-CRE**
BEDs (Supplementary Data 1–2), set `BACA_NECRE` / `BACA_ADCRE` — otherwise step 02
derives equivalents from the data.

## Run (in order)
```bash
source config.sh
python 00_organize.py            # builds work/manifest.tsv from the untarred files
#  >>> REVIEW work/manifest.tsv <<<  LuCaP NEPC/adeno is now auto-tagged (NE = 49/93/
#      145.1/145.2/173.1/208.1 -> fg; all other LuCaP -> PRAD/bg). Spot-check that split
#      and that `build` is right for each dataset (drives liftover); fix any condition=UNKNOWN.
bash   01_liftover.sh            # hg19 -> hg38 (BED via liftOver, bigWig via CrossMap)
export MANIFEST=$WORK/manifest.hg38.tsv
bash   02_reference_regions.sh   # NE-CRE vs Ad-CRE region sets (hg38)
bash   03_motif_enrichment.sh    # HOMER known-motif enrichment, NEPC vs adeno
bash   04_signal_profiles.sh     # deepTools heatmap/profile of ChIP signal over NE/Ad-CRE
python 05_make_s15b_chromatin.py --outdir $OUT/figures   # the S15B motif panel
```

## Outputs
- `out/motif/*/knownResults.txt` — per-contrast HOMER enrichment.
- `out/figures/s15b_chromatin_motif.{pdf,png}` + `_data.tsv` — **the panel**: regulon
  TFs (rows) × NEPC-vs-adeno contrasts (cols), color = −log10 p. ASCL1/NKX2-1/NEUROD1/
  FOXA2 lighting up here = chromatin-level corroboration of the SCENIC/decoupleR regulons.
- `out/signal/{heatmap,profile}_NEvsAd.pdf` — FOXA1/FOXA2/H3K27ac/ASCL1 signal enriched
  at NE-CRE vs Ad-CRE (the companion signal evidence).

## What this does and does NOT show
- **Does:** the regulon TFs' motifs are over-represented in independent NEPC chromatin,
  and NEPC TF ChIP signal concentrates at NE-specific regulatory elements — across 2–3
  independent cohorts/models, in clinical-adjacent (PDX/patient) material.
- **Does NOT:** prove accessibility *in LTL331* (no matched ATAC) or causal TF binding.
  This is orthogonal, chromatin-consistent support, not a substitute for matched
  accessibility profiling.

## Gotchas
- **Build / liftover:** never run step 03/04 on a mix of hg19 peaks and hg38 regions — step 01
  must run first. The `build` column in `00_organize.py`'s `DATASETS` decides what gets lifted;
  it's a hardcoded guess (Baca/Wang hg19, Formaggio hg38). **Confirm with `bigWigInfo <file>.bw`
  before lifting** — a wrong hg38→"lift" silently shifts every coordinate. Liftover skips only
  `use==ignore` rows, so all fg/bg/signal hg19 files (incl. the adeno bg peaks) are lifted.
- **Manifest condition tags:** LuCaP NEPC/adeno is now auto-assigned (NE = 49/93/145.1/145.2/
  173.1/208.1 → fg; every other LuCaP → PRAD → bg). Earlier this left adeno as UNKNOWN→ignore,
  so **bg was empty and the NEPC-vs-adeno contrast silently collapsed** — spot-check the split
  in `manifest.tsv` anyway; garbage in `condition`/`use` → garbage contrasts.
- **HOMER `-nomotif`** in step 03 = known-motifs only (fast, sufficient for the TF-motif test).
  Remove it to also get de-novo motifs (slower; nice for the supplement).
- **CrossMap CLI:** v0.6+ is `CrossMap bigwig …`; older is `CrossMap.py bigwig …`. Step 01
  auto-detects.
