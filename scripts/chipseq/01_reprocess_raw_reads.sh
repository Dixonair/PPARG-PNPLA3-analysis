#!/usr/bin/env bash
# ==============================================================================
# FULL REPROCESSING of GSE115421 PPARG ChIP-seq: SRA reads -> hg38 BAM -> peaks
# Run on a Linux/macOS machine or HPC.
#
# ROUGH RESOURCES (approximate, not measured): ~30 Gbases of sequence in total,
# 16 GB RAM is comfortable, 100+ GB free disk for FASTQ/BAM, hours not minutes.
#
# TOOLS (bioconda):
#   conda create -n chip -c conda-forge -c bioconda sra-tools fastp bowtie2 samtools bedtools macs2
#   conda activate chip
#
# GENOME BUILD: this script aligns to GRCh38 (hg38). The authors' GEO files (bigWigs,
# HOMER peaks) are hg19, so never overlay them on hg38 results without a liftOver.
# Choose ONE build for everything you compare; if you prefer to match the authors and
# GSE41578 (also hg19), use a Bowtie2 hg19 index, ENCODE's hg19 blacklist, and set
# `build` to "hg19" in chipseq_step1_peaks_to_genes.R.
#
# IMPORTANT LIMITATION: the SRA run table for GSE115421 lists ten PPARG ChIP runs
# and NO input/IgG control, so peaks are called WITHOUT a control. That is less
# specific than a controlled call; keep the 3-of-5 consensus in step 1 and say so
# in Methods.
# ==============================================================================
# ---- guard: this script must never be able to close your terminal ------------------
# `source script.sh` or `. script.sh` runs it INSIDE your own shell, so every `exit`
# (and every error under `set -e`) would close the terminal. Always start it with:
#     bash scripts/reprocess_GSE115421_chipseq.sh [N]
if [ -z "${BASH_VERSION:-}" ]; then
  echo "Run this with bash:  bash $0" >&2; return 1 2>/dev/null || exit 1
fi
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  echo "This script was SOURCED, which would close your terminal on any error." >&2
  echo "Run it instead with:  bash ${BASH_SOURCE[0]} [N]" >&2
  return 1
fi

set -euo pipefail

# Run this from the PROJECT ROOT (e.g. ~/RProjects/ChIPseq). Layout used:
#   data/chip_reprocess/  working files (fastq temp, bam, peaks, logs)
#   data/ref/             Bowtie2 index + blacklist (large, reused)
#   data/peaks/           final GSM*_peaks.narrowPeak copies (read by the R script)
ROOT="$(pwd)"
# Everything printed is also saved to a log, so if a session drops you can read what happened.
LOG="$ROOT/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR
echo "Log file: $LOG"
THREADS="${THREADS:-4}"            # shared server: be modest, e.g.  THREADS=8 bash ...
WORK="${WORK:-$ROOT/data/chip_reprocess}"
REF="${REF:-$ROOT/data/ref}"
PEAKS_OUT="${PEAKS_OUT:-$ROOT/data/peaks}"
mkdir -p "$WORK"/{fastq,bam,peaks,logs} "$REF" "$PEAKS_OUT"
cd "$WORK"

# ---- sample sheet: GSM, SRR, group  (taken from your SraRunTable.csv) -----------
cat > samples.tsv << 'EOF'
GSM3177950	SRR7271363	hASC_adipocytes
GSM3177951	SRR7271364	hASC_adipocytes
GSM3177952	SRR7271365	hASC_adipocytes
GSM3177953	SRR7271366	hASC_adipocytes
GSM3177954	SRR7271367	hASC_adipocytes
GSM3416379	SRR7963286	adipose_tissue
GSM3416380	SRR7963287	adipose_tissue
GSM3416381	SRR7963288	adipose_tissue
GSM3416382	SRR7963289	adipose_tissue
GSM3416383	SRR7963290	adipose_tissue
EOF
# Optional argument = how many samples to process (first N lines of the sheet above):
#   bash reprocess_GSE115421_chipseq.sh 1     # TEST: only GSM3177950 (do this first)
#   bash reprocess_GSE115421_chipseq.sh       # all ten
if [ "${1:-}" != "" ]; then head -n "$1" samples.tsv > samples.run.tsv; else cp samples.tsv samples.run.tsv; fi
echo "Processing $(wc -l < samples.run.tsv) sample(s)"

# ---- reference (download ONCE; verify the links on the pages named) --------------
# Bowtie2 GRCh38 no-alt index: see https://benlangmead.github.io/aws-indexes/bowtie
#   (file GRCh38_noalt_as.zip) -> unzip into <project>/data/ref/GRCh38_noalt_as/
# ENCODE hg38 blacklist v2 (this URL was checked: 636 regions)
if [ ! -f "$REF/hg38-blacklist.v2.bed.gz" ]; then
  curl -L -o "$REF/hg38-blacklist.v2.bed.gz" \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
fi

BT2_INDEX="$REF/GRCh38_noalt_as/GRCh38_noalt_as"
if [ ! -f "${BT2_INDEX}.1.bt2" ]; then
  echo "Bowtie2 index not found at ${BT2_INDEX}.*  -- download it first (see comment above)." >&2; exit 1
fi
# ---- per-sample processing -------------------------------------------------------
while IFS=$'\t' read -r GSM SRR GROUP; do
  echo "=== $GSM ($SRR, $GROUP) ==="
  if [ -s "peaks/${GSM}_peaks.narrowPeak" ]; then echo "already done, skipping"; continue; fi

  # 1. download reads (single-end libraries: one FASTQ per run)
  prefetch "$SRR" -O fastq --max-size 100G
  # Pass the BARE ACCESSION, not a constructed path. prefetch and fasterq-dump
  # share the same SRA toolkit cache/resolution logic (the config set up via
  # `vdb-config --interactive`); letting fasterq-dump resolve "$SRR" itself is
  # the standard, documented way to chain these two tools, and avoids having
  # to guess or reconstruct wherever prefetch actually cached the file (which
  # earlier attempts got wrong -- prefetch's "-O fastq" does not necessarily
  # mean the file ends up findable as a plain path under fastq/$SRR/).
  fasterq-dump "$SRR" -e "$THREADS" -O fastq --temp fastq
  gzip -f "fastq/${SRR}.fastq"

  # 2. adapter/quality trimming (report kept in logs/)
  fastp -i "fastq/${SRR}.fastq.gz" -o "fastq/${GSM}.trim.fastq.gz" -w "$THREADS" \
        -j "logs/${GSM}.fastp.json" -h "logs/${GSM}.fastp.html"

  # 3. align to GRCh38, keep MAPQ>=30, sort, index
  bowtie2 -p "$THREADS" -x "$BT2_INDEX" -U "fastq/${GSM}.trim.fastq.gz" 2> "logs/${GSM}.bowtie2.log" \
    | samtools view -@ "$THREADS" -b -q 30 - \
    | samtools sort -@ "$THREADS" -o "bam/${GSM}.q30.bam" -
  samtools index "bam/${GSM}.q30.bam"

  # 4. remove ENCODE blacklist regions
  bedtools intersect -v -abam "bam/${GSM}.q30.bam" -b "$REF/hg38-blacklist.v2.bed.gz" > "bam/${GSM}.filt.bam"
  samtools index "bam/${GSM}.filt.bam"

  # 5. peaks WITHOUT a control (none exists in the SRA table); MACS2 default
  #    duplicate handling. If model building fails, retry with a fixed fragment size.
  macs2 callpeak -t "bam/${GSM}.filt.bam" -f BAM -g hs -n "$GSM" -q 0.01 --outdir peaks \
        2> "logs/${GSM}.macs2.log" \
  || { echo "MACS2 model failed for $GSM; retrying with --nomodel --extsize 200" | tee -a "logs/${GSM}.macs2.log"
       macs2 callpeak -t "bam/${GSM}.filt.bam" -f BAM -g hs -n "$GSM" -q 0.01 --nomodel --extsize 200 \
             --outdir peaks 2>> "logs/${GSM}.macs2.log"; }

  cp "peaks/${GSM}_peaks.narrowPeak" "$PEAKS_OUT/"     # copy where the R script looks

  # 6. free disk: raw reads are large and can be re-downloaded
  rm -rf "fastq/$SRR" "fastq/${SRR}.fastq.gz" "fastq/${GSM}.trim.fastq.gz"

  echo "$GSM: $(wc -l < peaks/${GSM}_peaks.narrowPeak) peaks; unique-mapped summary in logs/${GSM}.bowtie2.log"
done < samples.run.tsv

echo
echo "Done. Peaks: $WORK/peaks/ (originals) and $PEAKS_OUT/ (copies)."
echo "Copies are already in $PEAKS_OUT (read by chipseq_step1_peaks_to_genes.R)."
echo "Sanity checks before trusting anything: (1) alignment rate in logs/*.bowtie2.log,"
echo "(2) peak counts per sample (tens of thousands is typical for PPARG in human adipocytes),"
echo "(3) look at PNPLA3 (chr22, hg38) and FABP4 in a genome browser."
