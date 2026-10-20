#!/usr/bin/env bash
# Step 4: pooled-replicate peak calling for HepG2 PPARG ChIP-seq (ENCODE ENCSR130VQL), plus relaxed (q<0.05)
# per-replicate calls, pooled tracks and a pooled QC line. Needs step 1 finished (replicate BAMs + control BAM).
# Run from the project root with the conda env active.
#
# Optional environment variables
#   OUT       results/hepg2_ppar
#   CTRL_BAM  $OUT/bam/control_sub20.bam
#   THREADS   8
#   MACS      macs3 (default; falls back to macs2)
#   FORCE=1   redo everything
set -euo pipefail
trap 'echo "FAILED at line $LINENO (see $OUT/qc/*.log)" >&2' ERR
OUT="${OUT:-results/hepg2_ppar}"; THREADS="${THREADS:-8}"
CTRL_BAM="${CTRL_BAM:-$OUT/bam/control_sub20.bam}"
MACS="${MACS:-$(command -v macs3 || command -v macs2)}"
FORCE="${FORCE:-0}"; GSIZE=2913022398
mkdir -p "$OUT"/{peaks/relaxed,tracks,qc}
[[ -s "$CTRL_BAM" ]] || { echo "Control BAM not found: $CTRL_BAM" >&2; exit 1; }

# replicate BAMs = the ChIP samples, i.e. those with a per-replicate peak file from step 1 (this excludes the input/control BAMs).
# Override with REPS="ENCFF764SMY,ENCFF660VHF" if needed.
BAMS=()
if [[ -n "${REPS:-}" ]]; then IFS=',' read -ra LABS <<< "$REPS"
else
  LABS=(); for pk in "$OUT"/peaks/*_peaks.narrowPeak; do n=$(basename "$pk" _peaks.narrowPeak); [[ "$n" == ChIP_pooled* ]] && continue; LABS+=("$n"); done
fi
for n in "${LABS[@]}"; do [[ -s "$OUT/bam/$n.bam" ]] && BAMS+=("$OUT/bam/$n.bam"); done
(( ${#BAMS[@]} >= 2 )) || { echo "Need >=2 replicate BAMs in $OUT/bam (found ${#BAMS[@]})" >&2; exit 1; }
echo "Replicate BAMs: ${BAMS[*]}"

POOL="$OUT/bam/ChIP_pooled.bam"
if [[ "$FORCE" == 1 || ! -s "$POOL" ]]; then
  samtools merge -f -@ "$THREADS" "$POOL" "${BAMS[@]}"; samtools index "$POOL"
fi

call () {   # $1=bam  $2=name  $3=qvalue  $4=outdir
  local out="$4/$2_peaks.narrowPeak"
  if [[ "$FORCE" != 1 && -s "$out" ]]; then echo "[skip] $2 q<$3" >&2; return; fi
  "$MACS" callpeak -t "$1" -c "$CTRL_BAM" -f BAM -g hs -q "$3" --keep-dup all -n "$2" --outdir "$4" 2> "$OUT/qc/$2.macs.log"
}

call "$POOL" ChIP_pooled 0.01 "$OUT/peaks"
call "$POOL" ChIP_pooled_q05 0.05 "$OUT/peaks/relaxed"
for b in "${BAMS[@]}"; do lab=$(basename "$b" .bam); call "$b" "${lab}_q05" 0.05 "$OUT/peaks/relaxed"; done

if [[ "$FORCE" == 1 || ! -s "$OUT/tracks/ChIP_pooled.RPGC.bw" ]]; then
  bamCoverage -b "$POOL" -o "$OUT/tracks/ChIP_pooled.RPGC.bw" --normalizeUsing RPGC --effectiveGenomeSize "$GSIZE" \
      --binSize 10 --extendReads 200 -p "$THREADS"
fi
if [[ "$FORCE" == 1 || ! -s "$OUT/tracks/ChIP_pooled.log2ChIPoverInput.bw" ]]; then
  bamCompare -b1 "$POOL" -b2 "$CTRL_BAM" --operation log2 --scaleFactorsMethod readCount --pseudocount 1 \
      --binSize 10 --extendReads 200 -p "$THREADS" -o "$OUT/tracks/ChIP_pooled.log2ChIPoverInput.bw"
fi
# the control track is needed by step 3/5
if [[ ! -s "$OUT/tracks/control.RPGC.bw" ]]; then
  bamCoverage -b "$CTRL_BAM" -o "$OUT/tracks/control.RPGC.bw" --normalizeUsing RPGC --effectiveGenomeSize "$GSIZE" \
      --binSize 10 --extendReads 200 -p "$THREADS"
fi

tot=$(samtools view -c "$POOL"); npk=$(wc -l < "$OUT/peaks/ChIP_pooled_peaks.narrowPeak")
inpk=$(samtools view -c -L "$OUT/peaks/ChIP_pooled_peaks.narrowPeak" "$POOL")
echo -e "ChIP_pooled\treads=$tot\tpeaks=$npk\tFRiP=$(python3 -c "print(round($inpk/$tot,4))")" | tee -a "$OUT/qc/summary.tsv"
echo "Done. Next: python3 scripts/chipseq_hepg2/05_hepg2_full_analysis.py"
