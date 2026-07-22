#!/usr/bin/env python3
"""
00_organize.py — inventory the untarred GEO files and build an EDITABLE manifest.

GEO _RAW.tar filenames are not standardized, so this does best-effort tagging from
the filename (factor, condition, type) and writes manifest.tsv. You then REVIEW it
(especially the `condition` and `use` columns) before running the rest of the
pipeline — the downstream scripts read this manifest, nothing else.

Columns:
  dataset    baca | wang | formaggio
  build      hg19 | hg38           (baca/wang = hg19, formaggio = hg38)
  type       peak | signal         (BED/broadPeak/narrowPeak vs bigWig)
  factor     FOXA1 | FOXA2 | AR | JUN | FOSL1 | H3K27ac | H3K4me2 | ASCL1 | ATAC | input | NA
  condition  NEPC | ARPC | DNPC | UNKNOWN | IGNORE          <-- CHECK THESE
  use        fg | bg | signal | dnpc_fg | arpc_bg | ignore  <-- CHECK THESE
  path       absolute path to the file

Conditions (subtypes, not just NE-vs-not):
  NEPC = AR-/NE+ (neuroendocrine)        ARPC = AR-driven adenocarcinoma
  DNPC = AR-/NE- double-negative (mesenchymal/basal; the human analog of cluster 17)

Use semantics:
  fg/bg      = NEPC(fg) vs adeno(bg) PEAKS pooled in step 02 to define NE-CRE/Ad-CRE
  signal     = bigWig to profile over those regions in step 04
  dnpc_fg/arpc_bg = the patient-matched Wang 201.2(DNPC)-vs-201.1(ARPC) H3K27ac peak
               contrast in step 03 (the cluster-17-bridge sub-analysis); kept SEPARATE
               from the NE-CRE pool so DNPC chromatin never blurs the terminal-NE signal
  ignore     = excluded (inputs, perturbation arms, unclassified)

Authoritative phenotype sources (no guessing):
  * Baca GSE161948 = the GEO series-matrix `histology:` field per LuCaP line.
  * Wang GSE232555 = the per-GSM sample titles (a mechanistic FOXA2/AP-1 study, so only
    BASELINE vehicle/siNC/parental tracks are allow-listed; every si*/OE/KO/drug arm is
    dropped). 201.1=ARPC(dura) & 201.2=DNPC(lung) are matched mets from one patient.

Run:  python 00_organize.py   (after editing + sourcing config.sh; reads env vars)
"""
import os, re, sys
from pathlib import Path

DATASETS = {
    'baca':      (os.environ.get('RAW_BACA', ''), 'hg19'),
    'wang':      (os.environ.get('RAW_WANG', ''), 'hg19'),
    'formaggio': (os.environ.get('RAW_FORM', ''), 'hg38'),
}
MANIFEST = os.environ.get('MANIFEST', 'manifest.tsv')

PEAK_EXT   = ('.bed', '.broadpeak', '.narrowpeak', '.bed.gz', '.broadpeak.gz', '.narrowpeak.gz')
SIGNAL_EXT = ('.bw', '.bigwig', '.bw.gz', '.bigwig.gz')

# --- Baca GSE161948: LuCaP lines that are Neuroendocrine per the GEO histology field.
#     (Verified against GSE161948_series_matrix `histology:`; every OTHER LuCaP = adeno.)
#     NB this CORRECTS the old regex guess: 173.2 is adeno (NOT NE), and 208.1 IS NE.
BACA_NE_LUCAP = {'49', '93', '145.1', '145.2', '173.1', '208.1'}

# --- Wang GSE232555: curated BASELINE allow-list (perturbation arms excluded by omission).
#     Maps GSM -> subtype; any Wang file whose GSM is not here -> condition IGNORE -> use ignore.
#     201.1=ARPC (dura met), 201.2=DNPC (lung met) — patient-matched; H660=NEPC; PC3=DNPC.
WANG_GSM = {
    # H660 (NEPC) baseline — vehicle / siNC
    'GSM7350173': 'NEPC',  # FOXA2 vehicle
    'GSM7350175': 'NEPC',  # JUN siNC
    'GSM7350176': 'NEPC',  # H3K4me2 vehicle
    'GSM7350177': 'NEPC',  # ATAC vehicle
    # PDX201.2 (DNPC) baseline
    'GSM7350201': 'DNPC', 'GSM7350202': 'DNPC', 'GSM7350203': 'DNPC', 'GSM7350204': 'DNPC',  # FOXA2
    'GSM7350205': 'DNPC', 'GSM7350206': 'DNPC',                                              # H3K27ac
    # PC3 (DNPC) baseline — vehicle / siNC
    'GSM7350184': 'DNPC',  # FOXA2 vehicle
    'GSM7350187': 'DNPC',  # FOSL1 siNC
    'GSM7350188': 'DNPC',  # JUN siNC
    'GSM7350189': 'DNPC',  # H3K27ac vehicle
    'GSM7350190': 'DNPC',  # ATAC vehicle
    # PDX201.1 (ARPC) baseline
    'GSM7350192': 'ARPC', 'GSM7350193': 'ARPC',                                              # AR
    'GSM7350194': 'ARPC', 'GSM7350195': 'ARPC', 'GSM7350196': 'ARPC', 'GSM7350197': 'ARPC',  # FOXA1
    'GSM7350198': 'ARPC', 'GSM7350199': 'ARPC', 'GSM7350200': 'ARPC',                        # H3K27ac
}

def tag_factor(name):
    n = name.upper()
    if 'FOXA1' in n: return 'FOXA1'
    if 'FOXA2' in n: return 'FOXA2'
    if re.search(r'H3K27AC|H3K27_AC', n): return 'H3K27ac'
    if 'H3K4ME1' in n: return 'H3K4me1'
    if 'H3K4ME2' in n: return 'H3K4me2'
    if 'H3K4ME3' in n: return 'H3K4me3'
    if 'H3K27ME3' in n: return 'H3K27me3'
    if 'ASCL1' in n: return 'ASCL1'
    if 'ATAC' in n: return 'ATAC'
    if re.search(r'FOSL1|FRA1', n): return 'FOSL1'
    if re.search(r'(?:CHIP[-_]?)?\bJUN\b|[_-]JUN(?![A-Z])', n): return 'JUN'
    if re.search(r'(?:CHIP[-_]?AR|[_-]AR)(?![A-Z0-9])', n): return 'AR'
    if re.search(r'INPUT|IGG', n): return 'input'
    return 'NA'

def tag_condition(ds, name):
    n = name.upper()
    if ds == 'wang':
        # Baseline allow-list only; bare-numeric PDX IDs + perturbation arms make
        # filename heuristics unsafe here, so classify strictly by GSM accession.
        m = re.search(r'GSM\d+', n)
        return WANG_GSM.get(m.group(0), 'IGNORE') if m else 'IGNORE'
    if ds == 'baca':
        # Authoritative GEO histology: NE set above, every other LuCaP = adeno (ARPC).
        m = re.search(r'LUCAP[_ ]?([0-9]+(?:\.[0-9])?)(?:CR)?', n)
        if m:
            return 'NEPC' if m.group(1) in BACA_NE_LUCAP else 'ARPC'
        return 'UNKNOWN'
    if ds == 'formaggio':
        # bigWig-only NEPC signal context (H660, LuCaP145 PDX, MSKPCa1).
        if re.search(r'H660|145|MSKPCA1|NEPC', n): return 'NEPC'
        return 'UNKNOWN'
    # generic fallback for any other dataset
    if re.search(r'NCIH660|NCI_H660|H660|42DENZR|NEPC', n): return 'NEPC'
    if re.search(r'LNCAP|22RV1|VCAP', n): return 'ARPC'
    if re.search(r'\bPC3\b|PC_3|DU145', n): return 'DNPC'
    return 'UNKNOWN'

def default_use(ds, ftype, factor, cond):
    if cond in ('IGNORE', 'UNKNOWN'):
        return 'ignore'
    if ftype == 'signal':
        prof = ('FOXA1', 'FOXA2', 'AR', 'JUN', 'FOSL1', 'H3K27ac', 'H3K4me2', 'ASCL1', 'ATAC')
        return 'signal' if factor in prof else 'ignore'
    if ftype == 'peak':
        if factor == 'input': return 'ignore'
        if ds == 'wang':
            # Only the patient-matched 201.2(DNPC)-vs-201.1(ARPC) H3K27ac pair becomes a
            # peak contrast (the cl17-bridge sub-analysis). H660/PC3 cell-line peaks stay
            # signal-only — never pooled into the Baca-anchored NE-CRE region set.
            if factor == 'H3K27ac' and cond == 'DNPC': return 'dnpc_fg'
            if factor == 'H3K27ac' and cond == 'ARPC': return 'arpc_bg'
            return 'ignore'
        if cond == 'NEPC': return 'fg'
        if cond == 'ARPC': return 'bg'
        return 'ignore'   # DNPC / other -> not part of the NE-CRE contrast
    return 'ignore'

def main():
    rows = []
    for ds, (root, build) in DATASETS.items():
        if not root or not Path(root).is_dir():
            print(f"[warn] {ds}: dir not found ({root!r}) — skipping", file=sys.stderr)
            continue
        for p in sorted(Path(root).rglob('*')):
            if not p.is_file():
                continue
            low = p.name.lower()
            if low.endswith(PEAK_EXT):
                ftype = 'peak'
            elif low.endswith(SIGNAL_EXT):
                ftype = 'signal'
            else:
                continue
            factor = tag_factor(p.name)
            cond   = tag_condition(ds, p.name)
            use    = default_use(ds, ftype, factor, cond)
            rows.append([ds, build, ftype, factor, cond, use, str(p)])

    Path(MANIFEST).parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, 'w') as fh:
        fh.write('dataset\tbuild\ttype\tfactor\tcondition\tuse\tpath\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')

    from collections import Counter
    uc = Counter(r[5] for r in rows)
    n_unk = sum(1 for r in rows if r[4] == 'UNKNOWN')
    print(f"wrote {MANIFEST}: {len(rows)} files  "
          f"(fg={uc['fg']}, bg={uc['bg']}, signal={uc['signal']}, "
          f"dnpc_fg={uc['dnpc_fg']}, arpc_bg={uc['arpc_bg']}, "
          f"ignore={uc['ignore']}, UNKNOWN-condition={n_unk})")
    print("ACTION: spot-check manifest.tsv before continuing —")
    print("  • Baca/Formaggio: any condition=UNKNOWN rows (unmatched LuCaP/sample naming)")
    print("  • Wang: confirm the allow-listed GSMs resolved (fg/bg should be 0 for wang; it")
    print("    contributes signal + the dnpc_fg/arpc_bg matched pair only)")

if __name__ == '__main__':
    main()
