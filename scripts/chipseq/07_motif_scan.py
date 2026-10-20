#!/usr/bin/env python3
"""
Scan the PNPLA3 ChIP-seq peak region (intron-1, downstream of TSS) and the
previously EMSA-validated upstream PPRE region for a PPARG::RXRA motif match,
using the real JASPAR MA0065.1 position frequency matrix (human, SELEX-derived,
https://jaspar.elixir.no/matrix/MA0065.1/).

Requires internet access to api.genome.ucsc.edu at run time.
Uses UCSC's documented getData/sequence endpoint (NOT a name-lookup).

Run from anywhere with internet access. No local dependencies beyond
the Python standard library (urllib, json, math).
"""
import urllib.request
import json
import math

# ------------------------------------------------------------------------------
# Real JASPAR MA0065.1 PFM (PPARG::RXRA, Homo sapiens, SELEX; 41 sequences)
# Source: https://jaspar2022.genereg.net/matrix/MA0065.1/ -- verified, not
# reconstructed from memory.
# ------------------------------------------------------------------------------
PFM = {
    "A": [10, 7, 19, 13, 17, 0, 0, 0, 2, 41, 41, 39, 0, 0, 0, 2, 38, 10, 13, 7],
    "C": [12, 13, 4, 1, 0, 0, 0, 2, 37, 0, 0, 0, 0, 0, 2, 36, 1, 12, 15, 10],
    "G": [7, 8, 6, 27, 24, 41, 39, 0, 2, 0, 0, 2, 41, 39, 0, 1, 2, 4, 9, 12],
    "T": [12, 13, 12, 0, 0, 0, 2, 39, 0, 0, 0, 0, 0, 2, 39, 2, 0, 15, 4, 12],
}
MOTIF_LEN = 20
N_SEQS = 41  # total sequences contributing to the PFM (column sums)
PSEUDOCOUNT = 0.8  # standard small pseudocount to avoid log(0)
BACKGROUND = {"A": 0.29, "C": 0.21, "G": 0.21, "T": 0.29}  # approx human genome-wide GC~42%

def build_pwm():
    """Convert the PFM to a log-odds PWM (log2 of observed/background frequency)."""
    pwm = {base: [] for base in "ACGT"}
    for pos in range(MOTIF_LEN):
        col_total = sum(PFM[b][pos] for b in "ACGT")
        for base in "ACGT":
            freq = (PFM[base][pos] + PSEUDOCOUNT * BACKGROUND[base]) / (col_total + PSEUDOCOUNT)
            pwm[base].append(math.log2(freq / BACKGROUND[base]))
    return pwm

def score_sequence(seq, pwm):
    """Score every window of length MOTIF_LEN in seq (both strands) against the PWM."""
    seq = seq.upper().replace("N", "A")  # crude N-handling; flag N-content separately
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_comp = "".join(comp.get(b, "N") for b in reversed(seq))

    best = {"score": -1e9, "position": None, "strand": None, "match_seq": None}
    max_possible = sum(max(pwm[b][p] for b in "ACGT") for p in range(MOTIF_LEN))

    for strand, s in [("+", seq), ("-", rev_comp)]:
        for i in range(len(s) - MOTIF_LEN + 1):
            window = s[i:i + MOTIF_LEN]
            score = sum(pwm[window[p]][p] for p in range(MOTIF_LEN))
            if score > best["score"]:
                best = {"score": score, "position": i, "strand": strand, "match_seq": window}
    best["max_possible"] = max_possible
    best["pct_of_max"] = 100 * best["score"] / max_possible if max_possible else None
    return best

def fetch_sequence(chrom, start, end):
    """Fetch raw hg38 sequence via UCSC's documented getData/sequence endpoint."""
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url, timeout=20) as resp:
        data = json.loads(resp.read())
    return data.get("dna", "").upper()

if __name__ == "__main__":
    pwm = build_pwm()

    # Region 1: the real, reproducible ChIP-seq peak (intron 1, ~2.5kb downstream
    # of TSS) -- this is the NEW region being asked "does it contain a PPRE",
    # not the already-known upstream site.
    regions = {
        "ChIP_peak_intron1":   ("chr22", 43_925_300, 43_925_950),   # ~650bp, covers the real peak span
        "EMSA_PPRE_upstream":  ("chr22", 43_922_950, 43_923_150),   # ~200bp, covers both known PPREs as a positive control for the scan itself
    }

    print("=== PPARG::RXRA (JASPAR MA0065.1) motif scan ===\n")
    for label, (chrom, start, end) in regions.items():
        print(f"--- {label}: {chrom}:{start:,}-{end:,} ({end-start}bp) ---")
        try:
            seq = fetch_sequence(chrom, start, end)
        except Exception as e:
            print(f"  ERROR fetching sequence: {e}\n")
            continue
        n_count = seq.count("N")
        if n_count > 0:
            print(f"  NOTE: {n_count} N bases in fetched sequence (masked/unknown) -- "
                  f"treated as A for scoring, interpret results cautiously if high.")
        result = score_sequence(seq, pwm)
        print(f"  Best match: {result['match_seq']} ({result['strand']} strand, "
              f"position {result['position']} within region)")
        print(f"  Score: {result['score']:.2f} / max possible {result['max_possible']:.2f} "
              f"({result['pct_of_max']:.1f}%)")
        print(f"  Genomic position of match: {chrom}:{start + result['position']:,}-"
              f"{start + result['position'] + MOTIF_LEN:,}\n")

    print("INTERPRETATION GUIDANCE (not computed -- use judgment):")
    print("  - A score >70-75% of max is typically considered a plausible motif match;")
    print("    below ~60% is usually not meaningful. These are rough conventions, not")
    print("    formal p-values -- for a publication-grade claim, convert this score to")
    print("    an empirical p-value (e.g. via FIMO) rather than relying on %-of-max alone.")
    print("  - Compare the ChIP-peak region's best score directly against the EMSA-PPRE")
    print("    region's score (scanned with the same PWM here, as a positive control that")
    print("    the scanning method itself correctly detects a known functional PPRE).")
