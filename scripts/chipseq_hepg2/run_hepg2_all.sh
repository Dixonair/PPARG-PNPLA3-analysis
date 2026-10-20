#!/usr/bin/env bash
# Runs the whole HepG2 reanalysis after step 1 (alignment + per-replicate peaks) has finished:
#   00 rebuild the RNA-seq input tables from the original files (GSE162276 xlsx, ASC_Rosi.csv, bound_genes_promoter.csv)
#   04 pooled peaks  ->  03 Fig.6-style figure  ->  05 full pipeline (QC, input tests, motif, Hallmark, integration)
#   ->  07 revalidation (MACS peak-level statistics, local background, ortholog-based enrichment, formal between-dataset tests)  ->  06 manuscript figures
# Optional: ORTHO_FILE, MAP_BW, BLACKLIST_BED, EXTRA_BW (see the header of scripts/chipseq_hepg2/07_hepg2_revalidation.py)
# Run from the project root ("ChIPseq HepG2 analysis") with the conda env active:
#     bash scripts/chipseq_hepg2/run_hepg2_all.sh
# It first checks every input and tool and stops with a clear message if something is missing.
set -uo pipefail
SCR="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-results/hepg2_ppar}"
find_file () { for d in inputs ref data . ; do [[ -s "$d/$1" ]] && { echo "$d/$1"; return; }; done; }
GTF="${GTF:-$(find_file gencode.v46.basic.annotation.gtf.gz)}"
FASTA="${FASTA:-$(find_file hg38.analysisSet.fa.gz)}"
# originals for step 00 (set SRC=<folder> if they are somewhere else)
SRC="${SRC:-}"
if [[ -z "$SRC" ]]; then for d in inputs ref data . ; do [[ -s "$d/GSE162276_HF_TZD_CvsHF_C_DEG.xlsx" && -s "$d/ASC_Rosi.csv" && -s "$d/bound_genes_promoter.csv" ]] && { SRC="$d"; break; }; done; fi
DEST="${DEST:-inputs}"
DE_CSV="${DE_CSV:-$DEST/GSE162276_TZDvsControl_DE_results.csv}"
HASC_CSV="${HASC_CSV:-$DEST/ASCRosi_DE_with_ChIP_annotation.csv}"
DE_RECOMP_CSV="${DE_RECOMP_CSV:-$DEST/GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv}"
HALLMARK_PKL="${HALLMARK_PKL:-$(find_file hallmark_sets.pkl)}"
CTRL_BAM="${CTRL_BAM:-$OUT/bam/control_sub20.bam}"
bad=0
chk () { if [[ -z "${2:-}" || ! -s "$2" ]]; then echo "  MISSING  $1  ($2)"; bad=1; else echo "  ok       $1  ($2)"; fi; }
echo "== inputs =="; chk "GTF" "$GTF"; chk "genome FASTA" "$FASTA"
if [[ -z "$SRC" ]]; then echo "  MISSING  originals: GSE162276_HF_TZD_CvsHF_C_DEG.xlsx, ASC_Rosi.csv, bound_genes_promoter.csv (put them in inputs/ or set SRC=<folder>)"; bad=1
else echo "  ok       original RNA-seq files in $SRC"; fi
chk "Hallmark gene sets" "$HALLMARK_PKL"; chk "control BAM" "$CTRL_BAM"
n=$(ls "$OUT"/peaks/*_peaks.narrowPeak 2>/dev/null | grep -vc 'ChIP_pooled'); echo "  ChIP replicates with peak files: $n (need 2)"; (( n >= 2 )) || bad=1
for r in $(ls "$OUT"/peaks/*_peaks.narrowPeak 2>/dev/null | grep -v pooled); do echo "  ok       peaks $r"; done
echo "== tools =="
for t in samtools bamCoverage bamCompare; do command -v $t >/dev/null && echo "  ok       $t" || { echo "  MISSING  $t"; bad=1; }; done
command -v macs3 >/dev/null || command -v macs2 >/dev/null && echo "  ok       macs" || { echo "  MISSING  macs3/macs2"; bad=1; }
python3 - <<'P' || bad=1
import importlib
for m in ("pyBigWig", "pysam", "pandas", "numpy", "scipy", "matplotlib", "statsmodels"):
    try: importlib.import_module(m); print("  ok      ", m)
    except Exception: print("  MISSING ", m, "-> pip install", "pybigwig" if m == "pyBigWig" else m)
    if m != "statsmodels":
        try: importlib.import_module(m)
        except Exception: raise SystemExit(1)
P
if (( bad )); then echo; echo "Fix the MISSING items above, then rerun. Nothing was run."; exit 1; fi
echo; echo "== step 00: build the RNA-seq tables from the original files =="
SRC="$SRC" DEST="$DEST" python3 "$SCR/00_build_inputs_from_originals.py" || { echo "step 00 failed"; exit 1; }
chk "liver TZD DE table (authors' values)" "$DE_CSV"; chk "hASC rosiglitazone table" "$HASC_CSV"; (( bad )) && exit 1
export OUT GTF FASTA DE_CSV HASC_CSV HALLMARK_PKL CTRL_BAM DE_RECOMP_CSV
echo; echo "== step 04: pooled peaks =="; bash "$SCR/04_hepg2_pooled_peaks.sh" || { echo "step 04 failed"; exit 1; }
echo; echo "== step 03: Fig. 6-style figure =="; python3 "$SCR/03_hepg2_figure_and_integration.py" || { echo "step 03 failed"; exit 1; }
echo; echo "== step 05: full analysis =="; python3 "$SCR/05_hepg2_full_analysis.py" || { echo "step 05 failed"; exit 1; }
echo; echo "== step 07: revalidation (peak-level stats, local background, ortholog-based enrichment, interaction tests) =="; python3 "$SCR/07_hepg2_revalidation.py" || { echo "step 07 failed (steps 03-05 outputs are intact; the 06 figures need 07 for panels d/i - see the message above)"; exit 1; }
if [[ "${RUN_SENSITIVITY:-1}" == "1" && -s "$DE_RECOMP_CSV" ]]; then
  echo; echo "== step 07 (sensitivity): same analysis with the DESeq2 re-run of the liver read counts =="
  DE_CSV="$DE_RECOMP_CSV" REVAL_TAG="_DESeq2rerun" NPERM="${NPERM_SENS:-2000}" python3 "$SCR/07_hepg2_revalidation.py" || echo "  sensitivity run failed (primary results are unaffected)"
fi
echo; echo "== step 06: manuscript figures =="; python3 "$SCR/06_hepg2_manuscript_figures.py" || { echo "step 06 failed"; exit 1; }
echo; echo "== step 08: Source Data workbook =="; python3 "$SCR/08_export_source_data.py" || echo "  step 08 failed (figures are unaffected)"
echo; echo "ALL DONE. Figures: $OUT/full_analysis/manuscript_figures/ (HepG2_main_figure.pdf, HepG2_supplementary_figure.pdf) and Source_Data_Fig4_SuppFig4_HepG2.xlsx. Send me:"
echo "  $OUT/hepg2_PPARG_figure.pdf   $OUT/hepg2_integration_and_QC.txt"
echo "  $OUT/full_analysis/hepg2_full_analysis_figure.pdf   $OUT/full_analysis/hepg2_full_analysis_report.txt"
echo "  $OUT/full_analysis/window_enrichment_tests.csv   $OUT/qc/summary.tsv"
echo "  $DEST/inputs_provenance.txt   $OUT/full_analysis/revalidation/revalidation_report.txt   $OUT/full_analysis/revalidation_DESeq2rerun/revalidation_report.txt"
