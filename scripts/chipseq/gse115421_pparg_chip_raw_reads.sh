#!/usr/bin/env bash
###############################################################################
# GSE115421 PPARG ChIP-seq, ALL TEN SAMPLES, ONE UNIFORM PIPELINE
# (raw reads -> hg38 BAM -> peaks -> signal tracks -> QC table -> per-group consensus)
#   5 x hASC-derived adipocytes   GSM3177950-54   (SRR7271363-67, 50-nt single-end)
#   5 x human adipose tissue      GSM3416379-83   (SRR7963286-90, 100-nt single-end)
# The sample sheet is taken from the SraRunTable.csv of the series (BioProject PRJNA474921). Every sample goes through
# exactly the same tools and parameters; nothing is set per group.
#
# This script consolidates the project's two earlier scripts - reprocess_GSE115421_chipseq.sh (SRA -> peaks) and
# 03_bam_to_bigwig.sh (BAM -> bigWig) - WITHOUT changing a command or parameter that affects a BAM, peak or track.
# Added only: tool-version log, a QC summary table (reads, alignment rate, peaks, FRiP), an optional read-length
# trim (OFF by default), the option to keep FASTQ files, and a peak-level consensus per group.
#
# STATUS: the alignment/peak/track commands are those of the earlier scripts, which were not run in the chat
# environment (no sequencing tools there). The authors must confirm that this is the version that produced the peak
# files used in the manuscript; if a different command or parameter was run, change it here before publishing.
# Run ONE sample first (bash gse115421_pparg_chip_raw_reads.sh 1) and look at the QC line before running all ten.
#
# KNOWN LIMITATIONS (state them in the Methods):
#  * the series has NO input/IgG control, so peaks are called WITHOUT a control (MACS2 without -c): less specific than
#    a controlled call. Downstream, chipseq_step1_peaks_to_genes.R uses a >= 3-of-5 rule per group (MIN_DONORS = 3).
#  * read length differs between the groups (50 vs 100 nt). Set TRIM_TO=50 to cut every read to 50 nt if a
#    length-matched analysis is wanted; by default no length trimming is done.
#  * the authors' GEO tracks are hg19; this pipeline aligns to GRCh38 and does not liftOver.
#  * this is NOT the same pipeline as the HepG2 ENCODE re-analysis (MACS3, bowtie2 --very-sensitive, duplicates
#    removed, input control). Uniform here means uniform across the ten GSE115421 samples.
#  * the bedtools consensus below is PEAK-level (>= 3 of 5 samples overlap). The gene-level >= 3-of-5 rule in
#    chipseq_step1_peaks_to_genes.R is a different, gene-based definition; do not describe them as the same thing.
#
# TOOLS (versions are written to $WORK/logs/versions.txt):
#   conda create -n chip -c conda-forge -c bioconda sra-tools fastp bowtie2 samtools bedtools "macs2<2.2.9" deeptools
#   (the project environment.yml pins macs2<2.2.9)
#
# USAGE (from the project root, with bash - never `source` it):
#   bash gse115421_pparg_chip_raw_reads.sh 1      # test on the first sample only
#   bash gse115421_pparg_chip_raw_reads.sh        # all ten
# Optional environment variables: THREADS (4)  WORK  REF  PEAKS_OUT  BW_DIR  TRIM_TO (unset)  KEEP_FASTQ=1  FORCE=1
###############################################################################
if [ -z "${BASH_VERSION:-}" ]; then echo "Run this with bash:  bash $0" >&2; return 1 2>/dev/null || exit 1; fi
if [ "${BASH_SOURCE[0]}" != "$0" ]; then echo "This script was SOURCED; run it with:  bash ${BASH_SOURCE[0]} [N]" >&2; return 1; fi
set -euo pipefail

ROOT="$(pwd)"
LOG="$ROOT/run_$(date +%Y%m%d_%H%M%S).log"; exec > >(tee -a "$LOG") 2>&1
trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR
echo "Log file: $LOG"
THREADS="${THREADS:-4}"; FORCE="${FORCE:-0}"; KEEP_FASTQ="${KEEP_FASTQ:-0}"; TRIM_TO="${TRIM_TO:-}"
WORK="${WORK:-$ROOT/data/chip_reprocess}"; REF="${REF:-$ROOT/data/ref}"
PEAKS_OUT="${PEAKS_OUT:-$ROOT/data/peaks}"; BW_DIR="${BW_DIR:-$ROOT/data/bigwig}"
GENOME_SIZE=2913022398                      # effective genome size, GRCh38 (deepTools)
for t in prefetch fasterq-dump fastp bowtie2 samtools bedtools macs2 bamCoverage python3; do
  command -v "$t" >/dev/null || { echo "missing tool: $t" >&2; exit 1; }
done
mkdir -p "$WORK"/{fastq,bam,peaks,logs} "$REF" "$PEAKS_OUT" "$BW_DIR"; cd "$WORK"

# ---------------------------------------------------------------- versions
{ echo "date: $(date -u +%FT%TZ)"
  echo "fasterq-dump: $(fasterq-dump --version 2>&1 | grep -i -m1 version || true)"
  echo "fastp: $(fastp --version 2>&1 | head -1)"; echo "bowtie2: $(bowtie2 --version | head -1)"
  echo "samtools: $(samtools --version | head -1)"; echo "bedtools: $(bedtools --version)"
  echo "macs2: $(macs2 --version 2>&1 | head -1)"; echo "deepTools: $(bamCoverage --version 2>&1 | head -1)"
  echo "TRIM_TO=${TRIM_TO:-none}  THREADS=$THREADS  GENOME_SIZE=$GENOME_SIZE"; } > logs/versions.txt

# ---------------------------------------------------------------- sample sheet: GSM, SRR, group (series SraRunTable.csv)
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
if [ "${1:-}" != "" ]; then head -n "$1" samples.tsv > samples.run.tsv; else cp samples.tsv samples.run.tsv; fi
echo "Processing $(wc -l < samples.run.tsv) sample(s)"

# ---------------------------------------------------------------- reference
# Bowtie2 GRCh38 no-alt analysis-set index (GRCh38_noalt_as.zip from https://benlangmead.github.io/aws-indexes/bowtie),
# unzipped into $REF/GRCh38_noalt_as/ ; the ENCODE hg38 blacklist v2 is downloaded once.
if [ ! -f "$REF/hg38-blacklist.v2.bed.gz" ]; then
  curl -L -o "$REF/hg38-blacklist.v2.bed.gz" https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
fi
BT2_INDEX="$REF/GRCh38_noalt_as/GRCh38_noalt_as"
[ -f "${BT2_INDEX}.1.bt2" ] || { echo "Bowtie2 index not found at ${BT2_INDEX}.*  -- download it first." >&2; exit 1; }

# ---------------------------------------------------------------- processing of ONE sample (identical for all ten)
process_sample () {
  local GSM="$1" SRR="$2"
  # 1. download reads (single-end: one FASTQ per run). Bare accession is passed to fasterq-dump (SRA toolkit cache).
  prefetch "$SRR" -O fastq --max-size 100G < /dev/null
  fasterq-dump "$SRR" -e "$THREADS" -O fastq --temp fastq < /dev/null
  gzip -f "fastq/${SRR}.fastq"

  # 2. adapter/quality trimming with fastp defaults (optional uniform read length: TRIM_TO)
  local LEN_OPT=(); if [ -n "$TRIM_TO" ]; then LEN_OPT=(--max_len1 "$TRIM_TO"); fi
  fastp -i "fastq/${SRR}.fastq.gz" -o "fastq/${GSM}.trim.fastq.gz" -w "$THREADS" ${LEN_OPT[@]+"${LEN_OPT[@]}"} \
        -j "logs/${GSM}.fastp.json" -h "logs/${GSM}.fastp.html" < /dev/null

  # 3. align to GRCh38 (bowtie2 default sensitivity), keep MAPQ >= 30, sort, index
  bowtie2 -p "$THREADS" -x "$BT2_INDEX" -U "fastq/${GSM}.trim.fastq.gz" 2> "logs/${GSM}.bowtie2.log" \
    | samtools view -@ "$THREADS" -b -q 30 - \
    | samtools sort -@ "$THREADS" -o "bam/${GSM}.q30.bam" -
  samtools index "bam/${GSM}.q30.bam"

  # 4. remove ENCODE hg38 blacklist regions
  bedtools intersect -v -abam "bam/${GSM}.q30.bam" -b "$REF/hg38-blacklist.v2.bed.gz" > "bam/${GSM}.filt.bam"
  samtools index "bam/${GSM}.filt.bam"

  # 5. peaks WITHOUT a control (none in the series); MACS2 default duplicate handling and fragment model.
  #    If model building fails, retry with a fixed fragment size (the retry is logged).
  macs2 callpeak -t "bam/${GSM}.filt.bam" -f BAM -g hs -n "$GSM" -q 0.01 --outdir peaks 2> "logs/${GSM}.macs2.log" \
  || { echo "MACS2 model failed for $GSM; retrying with --nomodel --extsize 200" | tee -a "logs/${GSM}.macs2.log"
       macs2 callpeak -t "bam/${GSM}.filt.bam" -f BAM -g hs -n "$GSM" -q 0.01 --nomodel --extsize 200 \
             --outdir peaks 2>> "logs/${GSM}.macs2.log"; }

  # 6. signal track: reads per genomic content, 10-bp bins, smoothed over 30 bp
  bamCoverage -b "bam/${GSM}.filt.bam" -o "$BW_DIR/${GSM}.bw" --normalizeUsing RPGC \
              --effectiveGenomeSize "$GENOME_SIZE" --binSize 10 --smoothLength 30 -p "$THREADS"

  # 7. free disk: reads can be re-downloaded (KEEP_FASTQ=1 keeps them)
  if [ "$KEEP_FASTQ" != 1 ]; then rm -rf "fastq/$SRR" "fastq/${SRR}.fastq.gz" "fastq/${GSM}.trim.fastq.gz"; fi
}

# ---------------------------------------------------------------- main loop (sample list on fd 3 so no tool can swallow it)
while IFS=$'\t' read -r -u 3 GSM SRR GROUP; do
  echo "=== $GSM ($SRR, $GROUP) ==="
  if [ "$FORCE" != 1 ] && [ -s "peaks/${GSM}_peaks.narrowPeak" ] && [ -s "$BW_DIR/${GSM}.bw" ]; then
    echo "already done, skipping (FORCE=1 to redo)"
  else
    process_sample "$GSM" "$SRR"
  fi
  cp "peaks/${GSM}_peaks.narrowPeak" "$PEAKS_OUT/"       # copy read by chipseq_step1_peaks_to_genes.R
done 3< samples.run.tsv

# ---------------------------------------------------------------- QC table (rebuilt from the logs/BAMs of every sample in this run)
printf "sample\tgroup\tsrr\tinput_reads\tafter_trim_reads\talignment_rate\tfiltered_reads\tpeaks_q0.01\tFRiP\n" > logs/qc_summary.tsv
while IFS=$'\t' read -r -u 3 GSM SRR GROUP; do
  read -r n_in n_trim < <(python3 -c "import json;d=json.load(open('logs/${GSM}.fastp.json'))['summary'];print(d['before_filtering']['total_reads'],d['after_filtering']['total_reads'])")
  ar=$(grep "overall alignment rate" "logs/${GSM}.bowtie2.log" | awk '{print $1}')
  tot=$(samtools view -c "bam/${GSM}.filt.bam")
  inpk=$(samtools view -c -L "peaks/${GSM}_peaks.narrowPeak" "bam/${GSM}.filt.bam")
  npk=$(wc -l < "peaks/${GSM}_peaks.narrowPeak")
  frip=$(python3 -c "print(round($inpk/$tot,4))")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$GSM" "$GROUP" "$SRR" "$n_in" "$n_trim" "$ar" "$tot" "$npk" "$frip" >> logs/qc_summary.tsv
done 3< samples.run.tsv
if command -v column >/dev/null; then column -t -s$'\t' logs/qc_summary.tsv; else cat logs/qc_summary.tsv; fi

# ---------------------------------------------------------------- peak-level consensus per group (>= 3 of 5), only if all 5 peak files exist
mkdir -p peaks/sorted
for G in hASC_adipocytes adipose_tissue; do
  FILES=(); ok=1
  for GSM in $(awk -F'\t' -v g="$G" '$3==g{print $1}' samples.tsv); do
    f="peaks/${GSM}_peaks.narrowPeak"
    if [ -s "$f" ]; then sort -k1,1 -k2,2n "$f" > "peaks/sorted/${GSM}.bed"; FILES+=("peaks/sorted/${GSM}.bed"); else ok=0; fi
  done
  if [ "$ok" = 1 ] && [ "${#FILES[@]}" -eq 5 ]; then
    bedtools multiinter -i "${FILES[@]}" | awk -F'\t' '$4>=3' | bedtools merge -i - > "peaks/consensus_${G}_ge3of5.bed"
    echo "consensus $G: $(wc -l < "peaks/consensus_${G}_ge3of5.bed") regions covered by >= 3 of 5 samples"
  else
    echo "consensus $G skipped (not all 5 peak files present)"
  fi
done

echo; echo "Done. Peaks: $WORK/peaks/ (copies in $PEAKS_OUT); tracks: $BW_DIR; QC: $WORK/logs/qc_summary.tsv; versions: $WORK/logs/versions.txt"
echo "Next: Rscript chipseq_step1_peaks_to_genes.R   (gene annotation and the >=3-of-5 rule per group)"
echo "Sanity checks: alignment rates, peaks and FRiP per sample (FRiP of a no-control call is not comparable to a controlled one); inspect PNPLA3 (chr22, hg38) and FABP4 in a browser."
