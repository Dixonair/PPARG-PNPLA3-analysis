#!/usr/bin/env bash
###############################################################################
# PPARG ChIP-seq in HepG2 (ENCODE ENCSR130VQL, CETCh-seq): RAW READS -> filtered BAMs -> peaks -> signal tracks.
#
# This single script replaces steps 01 and 04 of the project (same commands, same parameters, same output file names, so
# 02/03/05/06/07/08 run unchanged on its output).  Differences from 01 + 04: (i) the control is subsampled by an explicit,
# seeded command (see CONTROL below), (ii) tool versions, parameters and input md5 sums are written to $OUT/qc, (iii) extra
# log files (samtools markdup statistics, flagstat).  None of these changes alters a BAM, peak or track.
#
# INPUT   single-end FASTQ (.fastq or .fastq.gz).  Two ChIP replicates and one input (control) library.
# OUTPUT  $OUT/bam/<label>.bam                 filtered, de-duplicated replicate BAMs
#         $OUT/bam/ChIP_pooled.bam             merge of the replicate BAMs
#         $OUT/bam/control_sub<PCT>.bam        control BAM used for ALL peak calls and tracks (default control_sub20.bam)
#         $OUT/peaks/<label>_peaks.narrowPeak  per-replicate MACS3 peaks, q < 0.01 (also ChIP_pooled_peaks.narrowPeak)
#         $OUT/peaks/relaxed/*_q05_peaks.narrowPeak   q < 0.05 calls (each replicate and pooled)
#         $OUT/tracks/*.RPGC.bw, *.log2ChIPoverInput.bw, control.RPGC.bw
#         $OUT/qc/                             fastp, bowtie2, markdup, flagstat, MACS logs, versions.txt, run_parameters.txt, summary.tsv
#
# REQUIRED environment variables
#   INDEX  bowtie2 index prefix of the human reference (GRCh38/hg38), e.g. /data/ref/bt2/GRCh38
#   CHIP   comma-separated ChIP FASTQ files, one per replicate (>= 2)
#   CTRL   input/control: comma-separated FASTQ file(s), or ONE ready-made .bam (then no alignment is done and no
#          subsampling is applied - the BAM is used as given)
# OPTIONAL
#   OUT (results/hepg2_ppar)   THREADS (8)   FORCE=1 (redo everything)   MACS (macs3, falls back to macs2)
#   CTRL_FRACTION  fraction of the processed control reads to keep, 0 < f < 1  (default 0.20 -> control_sub20.bam).
#                  Set CTRL_FRACTION=1 to use the whole control.
#   CTRL_SEED      integer seed for the subsampling (default 42)
#
# EXAMPLE (replace the placeholders with the real ENCODE file accessions; check "Experiment"/"Controlled by" on each ENCODE
# file page - do not infer ChIP vs control from the order of the accessions):
#   INDEX=/data/ref/bt2/GRCh38 \
#   CHIP=data/<ChIP_rep1>.fastq.gz,data/<ChIP_rep2>.fastq.gz  CTRL=data/<input>.fastq.gz \
#   THREADS=8 OUT=results/hepg2_ppar  bash hepg2_ppar_chip_raw_reads.sh
#
# CONTROL
#   The control BAM used for the published figures was named control_sub20.bam (31,198,466 mapped reads).  The command that
#   produced it was not recorded.  This script makes the step explicit and reproducible:
#       samtools view -b -s <SEED>.<FRACTION> control.bam > control_sub20.bam
#   Because the original command is unknown, rerunning this script will give a control BAM with the same size class but not
#   necessarily the same reads; peak calls and window statistics can therefore differ slightly from the published numbers.
#   Report the control that was actually used (qc/run_parameters.txt records it).
###############################################################################
set -euo pipefail
trap 'echo "FAILED at line $LINENO (see $OUT/qc/*.log)" >&2' ERR

: "${INDEX:?set INDEX (bowtie2 index prefix)}"; : "${CHIP:?set CHIP (comma-separated ChIP FASTQs)}"; : "${CTRL:?set CTRL (control FASTQ(s) or one .bam)}"
THREADS="${THREADS:-8}"; OUT="${OUT:-results/hepg2_ppar}"; FORCE="${FORCE:-0}"
CTRL_FRACTION="${CTRL_FRACTION:-0.20}"; CTRL_SEED="${CTRL_SEED:-42}"
MACS="${MACS:-$(command -v macs3 || command -v macs2 || true)}"
GSIZE=2913022398            # GRCh38 effective genome size for deepTools RPGC (non-N bases, 2022 value for unique + multi-mapped)
for t in fastp bowtie2 samtools bamCoverage bamCompare python3; do command -v "$t" >/dev/null || { echo "missing tool: $t" >&2; exit 1; }; done
[[ -n "$MACS" ]] || { echo "missing tool: macs3 (or macs2)" >&2; exit 1; }
mkdir -p "$OUT"/{trim,bam,peaks/relaxed,tracks,qc}
IFS=',' read -ra CHIPS <<< "$CHIP"
(( ${#CHIPS[@]} >= 2 )) || { echo "CHIP must list >= 2 replicate FASTQs (pooled calls need replicates)" >&2; exit 1; }

# ---------------------------------------------------------------- record versions, parameters, input checksums
{
  echo "date: $(date -u +%FT%TZ)"; echo "host: $(uname -srm)"
  echo "fastp: $(fastp --version 2>&1 | head -1)"; echo "bowtie2: $(bowtie2 --version | head -1)"
  echo "samtools: $(samtools --version | head -1)"; echo "MACS: $("$MACS" --version 2>&1 | head -1)"
  echo "deepTools: $(bamCoverage --version 2>&1 | head -1)"; echo "python: $(python3 --version 2>&1)"
} > "$OUT/qc/versions.txt"
{
  echo "INDEX=$INDEX"; echo "CHIP=$CHIP"; echo "CTRL=$CTRL"; echo "THREADS=$THREADS"; echo "CTRL_FRACTION=$CTRL_FRACTION"; echo "CTRL_SEED=$CTRL_SEED"
  echo "GSIZE(deepTools)=$GSIZE"; echo "MACS_gsize=hs"; echo "reads_filter=MAPQ>=30, samtools -F 1804, duplicates removed (markdup -r)"
  echo "--- md5 of input files"
  for f in "${CHIPS[@]}" ${CTRL//,/ }; do [[ -f "$f" ]] && md5sum "$f"; done
} > "$OUT/qc/run_parameters.txt"

# ---------------------------------------------------------------- 1. trim, align, filter, de-duplicate
align_one () {                 # $1 = FASTQ  ->  prints the label; writes $OUT/bam/<label>.bam (+ .bai)
  local fq="$1" lab; lab=$(basename "$fq"); lab=${lab%%.f*}
  if [[ "$FORCE" != 1 && -s "$OUT/bam/$lab.bam" ]]; then echo "[skip] $lab" >&2; echo "$lab"; return; fi
  fastp -i "$fq" -o "$OUT/trim/$lab.fq.gz" -w "$THREADS" -h "$OUT/qc/$lab.fastp.html" -j "$OUT/qc/$lab.fastp.json" >&2   # fastp defaults: adapter auto-detection, quality filter
  bowtie2 --very-sensitive -p "$THREADS" -x "$INDEX" -U "$OUT/trim/$lab.fq.gz" 2> "$OUT/qc/$lab.bowtie2.log" \
    | samtools view -@ "$THREADS" -b -q 30 -F 1804 - \
    | samtools sort -@ "$THREADS" -o "$OUT/bam/$lab.sorted.bam" - >&2
  samtools markdup -r -s -@ "$THREADS" "$OUT/bam/$lab.sorted.bam" "$OUT/bam/$lab.bam" 2> "$OUT/qc/$lab.markdup.txt"
  samtools index "$OUT/bam/$lab.bam"; samtools flagstat "$OUT/bam/$lab.bam" > "$OUT/qc/$lab.flagstat.txt"
  rm -f "$OUT/bam/$lab.sorted.bam"; echo "$lab"
}
CHIP_LABS=(); for f in "${CHIPS[@]}"; do CHIP_LABS+=("$(align_one "$f")"); done

# ---------------------------------------------------------------- 2. control: align (if FASTQ), then seeded subsample
PCT=$(awk -v f="$CTRL_FRACTION" 'BEGIN{printf "%d", f*100 + 0.5}')
CTRL_BAM="$OUT/bam/control_sub${PCT}.bam"
if [[ "$CTRL" == *.bam ]]; then
  CTRL_BAM="$CTRL"; echo "[control] using the supplied BAM as given: $CTRL_BAM" | tee -a "$OUT/qc/run_parameters.txt"
elif [[ "$FORCE" == 1 || ! -s "$CTRL_BAM" ]]; then
  IFS=',' read -ra CTRLS <<< "$CTRL"; CTRL_LABS=(); for f in "${CTRLS[@]}"; do CTRL_LABS+=("$(align_one "$f")"); done
  if (( ${#CTRL_LABS[@]} > 1 )); then
    samtools merge -f -@ "$THREADS" "$OUT/bam/control_full.bam" $(for l in "${CTRL_LABS[@]}"; do echo "$OUT/bam/$l.bam"; done)
    CFULL="$OUT/bam/control_full.bam"
  else CFULL="$OUT/bam/${CTRL_LABS[0]}.bam"; fi
  samtools index "$CFULL"
  if awk -v f="$CTRL_FRACTION" 'BEGIN{exit !(f>=1)}'; then
    cp "$CFULL" "$CTRL_BAM"                       # CTRL_FRACTION=1: whole control
  else
    SARG=$(awk -v s="$CTRL_SEED" -v f="$CTRL_FRACTION" 'BEGIN{printf "%d%s", s, substr(sprintf("%.6f", f), 2)}')   # e.g. 42 + 0.20 -> 42.200000
    samtools view -b -@ "$THREADS" -s "$SARG" "$CFULL" > "$CTRL_BAM"   # samtools: integer part = seed, fractional part = fraction kept
  fi
  samtools index "$CTRL_BAM"
  { echo "control_full_reads=$(samtools view -c "$CFULL")"; echo "control_used_reads=$(samtools view -c "$CTRL_BAM")"
    echo "control_subsample_arg=-s ${SARG:-none}"; } | tee -a "$OUT/qc/run_parameters.txt"
fi
[[ -s "$CTRL_BAM" ]] || { echo "control BAM not found: $CTRL_BAM" >&2; exit 1; }
[[ -s "$CTRL_BAM.bai" || -s "${CTRL_BAM%.bam}.bai" ]] || samtools index "$CTRL_BAM"

# ---------------------------------------------------------------- 3. pooled replicates
POOL="$OUT/bam/ChIP_pooled.bam"
if [[ "$FORCE" == 1 || ! -s "$POOL" ]]; then
  samtools merge -f -@ "$THREADS" "$POOL" $(for l in "${CHIP_LABS[@]}"; do echo "$OUT/bam/$l.bam"; done); samtools index "$POOL"
fi

# ---------------------------------------------------------------- 4. peak calling (MACS3, narrow peaks, default model, input control)
call () {                      # $1 BAM  $2 name  $3 q-value  $4 outdir
  local out="$4/$2_peaks.narrowPeak"
  if [[ "$FORCE" != 1 && -s "$out" ]]; then echo "[skip] $2 q<$3" >&2; return; fi
  "$MACS" callpeak -t "$1" -c "$CTRL_BAM" -f BAM -g hs -q "$3" --keep-dup all -n "$2" --outdir "$4" 2> "$OUT/qc/$2.macs.log"
}
for l in "${CHIP_LABS[@]}"; do
  call "$OUT/bam/$l.bam" "$l" 0.01 "$OUT/peaks"                 # primary: q < 0.01, per replicate
  call "$OUT/bam/$l.bam" "${l}_q05" 0.05 "$OUT/peaks/relaxed"   # relaxed: q < 0.05, per replicate
done
call "$POOL" ChIP_pooled 0.01 "$OUT/peaks"                      # pooled replicates, q < 0.01
call "$POOL" ChIP_pooled_q05 0.05 "$OUT/peaks/relaxed"          # pooled replicates, q < 0.05

# ---------------------------------------------------------------- 5. signal tracks (10-bp bins, reads extended to 200 bp)
track () {                     # $1 BAM  $2 name
  [[ "$FORCE" != 1 && -s "$OUT/tracks/$2.RPGC.bw" ]] || bamCoverage -b "$1" -o "$OUT/tracks/$2.RPGC.bw" --normalizeUsing RPGC \
      --effectiveGenomeSize "$GSIZE" --binSize 10 --extendReads 200 -p "$THREADS"
  [[ "$2" == control || ( "$FORCE" != 1 && -s "$OUT/tracks/$2.log2ChIPoverInput.bw" ) ]] || bamCompare -b1 "$1" -b2 "$CTRL_BAM" --operation log2 \
      --scaleFactorsMethod readCount --pseudocount 1 --binSize 10 --extendReads 200 -p "$THREADS" -o "$OUT/tracks/$2.log2ChIPoverInput.bw"
}
for l in "${CHIP_LABS[@]}"; do track "$OUT/bam/$l.bam" "$l"; done
track "$POOL" ChIP_pooled; track "$CTRL_BAM" control

# ---------------------------------------------------------------- 6. summary: library size, peaks, FRiP
: > "$OUT/qc/summary.tsv"
for l in "${CHIP_LABS[@]}" ChIP_pooled; do
  bam="$OUT/bam/$l.bam"; pk="$OUT/peaks/${l}_peaks.narrowPeak"
  tot=$(samtools view -c "$bam"); inpk=$(samtools view -c -L "$pk" "$bam"); npk=$(wc -l < "$pk")
  echo -e "$l\treads=$tot\tpeaks_q0.01=$npk\tFRiP=$(python3 -c "print(round($inpk/$tot,4))")" | tee -a "$OUT/qc/summary.tsv"
done
echo -e "control\treads=$(samtools view -c "$CTRL_BAM")" | tee -a "$OUT/qc/summary.tsv"
echo "Done. Next: python3 02_hepg2_quantify_pnpla3.py ; 05_hepg2_full_analysis.py (see run_hepg2_all.sh)"
