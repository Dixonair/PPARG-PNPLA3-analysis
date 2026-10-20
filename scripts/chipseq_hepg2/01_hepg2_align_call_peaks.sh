#!/usr/bin/env bash
# PPARG ChIP-seq in HepG2 (ENCODE ENCSR130VQL): align, call peaks, make tracks.
# Run in the RStudio *Terminal* tab (not the R console).
#
# Required environment variables:
#   INDEX   bowtie2 index prefix for hg38/GRCh38 (e.g. /data/ref/bt2/GRCh38)
#   CHIP    comma-separated ChIP FASTQ file(s)  (the PPARG ChIP, one per replicate)
# Optional:
#   CTRL    comma-separated input/control FASTQ file(s), OR a single ready-made control .bam. If empty, MACS2 runs without control.
#   THREADS (default 8)   OUT (default results/hepg2_ppar)
#
# NOTE: superseded by hepg2_ppar_chip_raw_reads.sh (steps 01 + 04 + seeded control subsampling + version log). Kept for the record.
# Example (placeholders - use the real ENCODE accessions of the two ChIP replicates and of the input):
#   INDEX=/data/ref/bt2/GRCh38 CHIP=data/<ChIP_rep1>.fastq,data/<ChIP_rep2>.fastq CTRL=data/<input>.fastq bash 01_hepg2_align_call_peaks.sh
# !! Check on ENCODE which file is the PPARG ChIP and which is the control
# !! (file page -> "Experiment"/"Controlled by"). Do NOT infer it from the order of the accessions.
set -euo pipefail
trap 'echo "FAILED at line $LINENO (see $OUT/qc/*.log)" >&2' ERR
: "${INDEX:?set INDEX}"; : "${CHIP:?set CHIP}"
THREADS="${THREADS:-8}"; OUT="${OUT:-results/hepg2_ppar}"; CTRL="${CTRL:-}"
MACS="${MACS:-$(command -v macs3 || command -v macs2)}"   # MACS3 preferred (MACS2 can fail on newer glibc)
GSIZE=2913022398                       # GRCh38 effective genome size (deepTools)
mkdir -p "$OUT"/{trim,bam,peaks,tracks,qc}

align_one () {                         # $1 = fastq   -> $OUT/bam/<label>.bam
  local fq="$1" lab; lab=$(basename "$fq"); lab=${lab%%.f*}
  [[ -s "$OUT/bam/$lab.bam" ]] && { echo "[skip] $lab" >&2; echo "$lab"; return; }
  fastp -i "$fq" -o "$OUT/trim/$lab.fq.gz" -w "$THREADS" \
        -h "$OUT/qc/$lab.fastp.html" -j "$OUT/qc/$lab.fastp.json" >&2
  bowtie2 --very-sensitive -p "$THREADS" -x "$INDEX" -U "$OUT/trim/$lab.fq.gz" \
        2> "$OUT/qc/$lab.bowtie2.log" \
   | samtools view -@ "$THREADS" -b -q 30 -F 1804 - \
   | samtools sort -@ "$THREADS" -o "$OUT/bam/$lab.sorted.bam" - >&2
  samtools markdup -r -@ "$THREADS" "$OUT/bam/$lab.sorted.bam" "$OUT/bam/$lab.bam" >&2
  samtools index "$OUT/bam/$lab.bam" >&2
  rm -f "$OUT/bam/$lab.sorted.bam"
  echo "$lab"
}

IFS=',' read -ra CHIPS <<< "$CHIP"
CHIP_LABS=(); for f in "${CHIPS[@]}"; do CHIP_LABS+=("$(align_one "$f")"); done
CTRL_BAM=""
if [[ "$CTRL" == *.bam ]]; then
  CTRL_BAM="$CTRL"                      # a ready-made control BAM (e.g. subsampled)
elif [[ -n "$CTRL" ]]; then
  IFS=',' read -ra CTRLS <<< "$CTRL"; CTRL_LABS=()
  for f in "${CTRLS[@]}"; do CTRL_LABS+=("$(align_one "$f")"); done
  if (( ${#CTRL_LABS[@]} > 1 )); then
    samtools merge -f -@ "$THREADS" "$OUT/bam/control_pooled.bam" $(for l in "${CTRL_LABS[@]}"; do echo "$OUT/bam/$l.bam"; done)
    samtools index "$OUT/bam/control_pooled.bam"; CTRL_BAM="$OUT/bam/control_pooled.bam"
  else CTRL_BAM="$OUT/bam/${CTRL_LABS[0]}.bam"; fi
fi

for lab in "${CHIP_LABS[@]}"; do
  bam="$OUT/bam/$lab.bam"
  if [[ -n "$CTRL_BAM" ]]; then
    "$MACS" callpeak -t "$bam" -c "$CTRL_BAM" -f BAM -g hs -q 0.01 --keep-dup all -n "$lab" --outdir "$OUT/peaks" 2> "$OUT/qc/$lab.macs2.log"
    bamCompare -b1 "$bam" -b2 "$CTRL_BAM" --operation log2 --scaleFactorsMethod readCount \
        --pseudocount 1 --binSize 10 --extendReads 200 -p "$THREADS" -o "$OUT/tracks/$lab.log2ChIPoverInput.bw"
  else
    "$MACS" callpeak -t "$bam" -f BAM -g hs -q 0.01 --keep-dup all -n "$lab" --outdir "$OUT/peaks" 2> "$OUT/qc/$lab.macs2.log"
  fi
  bamCoverage -b "$bam" -o "$OUT/tracks/$lab.RPGC.bw" --normalizeUsing RPGC \
      --effectiveGenomeSize "$GSIZE" --binSize 10 --extendReads 200 -p "$THREADS"
  tot=$(samtools view -c "$bam"); inpk=$(samtools view -c -L "$OUT/peaks/${lab}_peaks.narrowPeak" "$bam")
  npk=$(wc -l < "$OUT/peaks/${lab}_peaks.narrowPeak")
  echo -e "$lab\treads=$tot\tpeaks=$npk\tFRiP=$(python3 -c "print(round($inpk/$tot,4))")" | tee -a "$OUT/qc/summary.tsv"
done
echo "Done. Next: python3 02_hepg2_quantify_pnpla3.py"
