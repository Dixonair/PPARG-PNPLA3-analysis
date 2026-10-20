#!/usr/bin/env bash
# ==============================================================================
# BAM -> bigWig -> signal-depth visualization inputs
# Run on YOUR server, inside the chip_pnpla3 conda environment.
#
# Run step by step and check each output before continuing.
# ==============================================================================
set -euo pipefail

if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  echo "This script was SOURCED, which could close your terminal on error." >&2
  echo "Run it instead with:  bash ${BASH_SOURCE[0]}" >&2
  return 1
fi

ROOT="$(pwd)"
BAM_DIR="$ROOT/data/chip_reprocess/bam"
BW_DIR="$ROOT/data/bigwig"
SIGNAL_DIR="$ROOT/results/signal_extraction"
mkdir -p "$BW_DIR" "$SIGNAL_DIR"

SAMPLES="GSM3177950 GSM3177951 GSM3177952 GSM3177953 GSM3177954 GSM3416379 GSM3416380 GSM3416381 GSM3416382 GSM3416383"

echo "=== Step 0: confirm deeptools is available ==="
which bamCoverage computeMatrix plotProfile 2>&1 || {
  echo "deeptools not found -- installing into this conda environment..." >&2
  conda install -c bioconda -c conda-forge deeptools -y
}

echo
echo "=== Step 1: BAM -> bigWig (RPGC-normalised, standard for ChIP-seq) ==="
# RPGC (reads per genomic content) normalises for sequencing depth, the
# standard choice for comparing signal across samples with different totals
# (you have real depth differences -- see run_all*.log read counts).
GENOME_SIZE=2913022398   # effective genome size, GRCh38 (standard deeptools value)
for GSM in $SAMPLES; do
  BAM="$BAM_DIR/${GSM}.filt.bam"
  OUT="$BW_DIR/${GSM}.bw"
  if [ -f "$OUT" ]; then echo "$GSM: bigWig already exists, skipping"; continue; fi
  if [ ! -f "$BAM" ]; then echo "MISSING BAM for $GSM at $BAM -- skipping" >&2; continue; fi
  echo "Converting $GSM..."
  bamCoverage -b "$BAM" -o "$OUT" \
    --normalizeUsing RPGC --effectiveGenomeSize "$GENOME_SIZE" \
    --binSize 10 --smoothLength 30 -p 4
done

echo
echo "=== Step 2: genome-wide average signal around ALL TSS (classic ChIP-seq QC/result plot) ==="
# Needs a TSS BED file. If you don't have one, this step will tell you so;
# a simple way to get one is via UCSC Table Browser (RefSeq TSS, hg38, BED).
TSS_BED="$ROOT/data/ref/hg38_tss.bed"
if [ ! -f "$TSS_BED" ]; then
  echo "No TSS BED file found at $TSS_BED -- skipping the genome-wide profile plot." >&2
  echo "(Locus-specific PNPLA3/FABP4 extraction below does not need this file.)" >&2
else
  BW_LIST=$(for GSM in $SAMPLES; do echo -n "$BW_DIR/${GSM}.bw "; done)
  computeMatrix reference-point -S $BW_LIST -R "$TSS_BED" \
    --referencePoint TSS -b 3000 -a 3000 --binSize 50 -p 4 \
    -o "$SIGNAL_DIR/matrix_all_TSS.gz"
  plotProfile -m "$SIGNAL_DIR/matrix_all_TSS.gz" \
    -o "$SIGNAL_DIR/PPARG_signal_around_all_TSS.pdf" \
    --perGroup --plotTitle "PPARG ChIP-seq signal around all TSS"
  echo "Genome-wide TSS profile plot written: $SIGNAL_DIR/PPARG_signal_around_all_TSS.pdf"
fi

echo
echo "=== Step 3: extract raw signal values at PNPLA3 and FABP4 (for locus-specific plots) ==="
python3 scripts/extract_signal.py

echo
echo "DONE. Upload these two small files for plotting:"
echo "  $SIGNAL_DIR/PNPLA3_signal_profile.csv"
echo "  $SIGNAL_DIR/FABP4_signal_profile.csv"
echo "(and PPARG_signal_around_all_TSS.pdf from step 2, if it was produced, to view directly)"
