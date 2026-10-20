#!/usr/bin/env python3
"""
Multi-gene PPARG ChIP-seq signal heatmap: PNPLA3, FABP4, ADIPOQ, CD36, PLIN1,
plus true negative-control genes (not expressed/open in adipocytes).

Requires internet access to api.genome.ucsc.edu at run time.
The UCSC query logic uses the properly documented region-based API
(api.genome.ucsc.edu/getData/track?genome=...;track=...;chrom=...;start=...;end=...),
*** NOT a name-lookup (UCSC's API does not support querying by gene name directly).
*** A rough genomic window (sourced from GeneCards/OMIM/NCBI) is queried, then the
*** actual matching transcript's real strand/TSS is read from UCSC's live response --
*** this makes the script robust to any imprecision in the rough window itself.

Run from the project root (same directory as data/bigwig/).
Requires: pyBigWig
"""
import pyBigWig
import numpy as np
import csv
import os
import urllib.request
import json

bw_dir = os.path.join(os.getcwd(), "data/bigwig")
out_dir = os.path.join(os.getcwd(), "results/signal_extraction")
os.makedirs(out_dir, exist_ok=True)

samples = ["GSM3177950", "GSM3177951", "GSM3177952", "GSM3177953", "GSM3177954",
           "GSM3416379", "GSM3416380", "GSM3416381", "GSM3416382", "GSM3416383"]

# Rough genomic windows (from GeneCards/OMIM/NCBI, hg38) -- only used to query the
# right region; the ACTUAL strand/TSS used for extraction comes from UCSC's live
# ncbiRefSeqCurated track data returned for that region, not from these numbers.
rough_windows = {
    "ADIPOQ": ("chr3",  186_838_000, 186_862_000),   # established PPARG target
    "CD36":   ("chr7",   80_365_000,  80_685_000),   # established PPARG target
    "PLIN1":  ("chr15",  89_660_000,  89_685_000),   # established PPARG target
    "MYOD1":  ("chr11",  17_715_000,  17_726_000),   # skeletal muscle TF -- closed in adipocytes
    "RHO":    ("chr3",  129_525_000, 129_540_000),   # rod photoreceptor -- closed in adipocytes
    "INS":    ("chr11",   2_158_000,   2_162_000),   # pancreatic beta-cell -- closed in adipocytes
}

def fetch_gene_tss(gene_symbol, chrom, start, end):
    """Query UCSC for ncbiRefSeqCurated records in this region, find the record
    whose name2 (gene symbol) matches, and return its real strand/TSS."""
    url = (f"https://api.genome.ucsc.edu/getData/track?genome=hg38;"
           f"track=ncbiRefSeqCurated;chrom={chrom};start={start};end={end}")
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        print(f"  ERROR querying UCSC for {gene_symbol}: {e}")
        return None

    records = data.get("ncbiRefSeqCurated", [])
    matches = [r for r in records if r.get("name2", "").upper() == gene_symbol.upper()]
    if not matches:
        print(f"  WARNING: no record with name2=={gene_symbol} found in {chrom}:{start}-{end}")
        print(f"  Gene symbols actually found in this window: "
              f"{sorted(set(r.get('name2','?') for r in records))}")
        return None
    # if multiple transcripts, take the longest (most likely canonical/MANE-like)
    rec = max(matches, key=lambda r: r["txEnd"] - r["txStart"])
    strand = rec["strand"]
    tss = rec["txStart"] if strand == "+" else rec["txEnd"]
    return {"chrom": chrom, "tss": int(tss), "strand": strand, "transcript": rec.get("name", "?")}

print("=== Resolving real TSS/strand from live UCSC data (not hardcoded) ===")
regions = {}
for gene, (chrom, start, end) in rough_windows.items():
    info = fetch_gene_tss(gene, chrom, start, end)
    if info is None:
        print(f"SKIPPING {gene} -- could not confirm coordinates from UCSC. "
              f"Check manually at https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{start}-{end}")
        continue
    print(f"  {gene}: {info['transcript']} -> {info['chrom']}:{info['tss']:,} "
          f"({info['strand']} strand)")
    regions[gene] = (info["chrom"], info["tss"] - 2000, info["tss"] + 2000)

if not regions:
    raise SystemExit("No gene coordinates confirmed -- cannot proceed. Check network access "
                      "to api.genome.ucsc.edu from this server.")

# ------------------------------------------------------------------------------
# Extract bigWig signal for each confirmed gene, each sample
# ------------------------------------------------------------------------------
print("\n=== Extracting bigWig signal ===")
n_bins = 200
all_matrices = {}

for gene, (chrom, start, end) in regions.items():
    all_values = {}
    for s in samples:
        bw_path = os.path.join(bw_dir, s + ".bw")
        if not os.path.exists(bw_path):
            print(f"MISSING bigWig for {s}, filling with NA for {gene}")
            all_values[s] = [None] * n_bins
            continue
        bw = pyBigWig.open(bw_path)
        if chrom not in bw.chroms():
            print(f"WARNING: {chrom} not found in {s}.bw -- filling 0 for {gene}")
            all_values[s] = [0] * n_bins
            bw.close()
            continue
        vals = bw.stats(chrom, start, end, nBins=n_bins, type="mean")
        all_values[s] = [v if v is not None else 0 for v in vals]
        bw.close()
    all_matrices[gene] = all_values
    print(f"{gene}: extracted {chrom}:{start}-{end}")

for gene, (chrom, start, end) in regions.items():
    bin_size = (end - start) / n_bins
    positions = [int(start + i * bin_size) for i in range(n_bins)]
    rows = [["position"] + samples]
    for i, pos in enumerate(positions):
        rows.append([pos] + [all_matrices[gene][s][i] for s in samples])
    out_csv = os.path.join(out_dir, f"{gene}_signal_profile.csv")
    with open(out_csv, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"Wrote {out_csv}")

summary_rows = [["gene", "sample", "group", "mean_signal"]]
for gene in regions:
    for s in samples:
        vals = [v for v in all_matrices[gene][s] if v is not None]
        mean_val = float(np.mean(vals)) if vals else None
        group = "hASC" if s.startswith("GSM3177") else "tissue"
        summary_rows.append([gene, s, group, mean_val])

summary_path = os.path.join(out_dir, "multigene_heatmap_summary.csv")
with open(summary_path, "w", newline="") as f:
    csv.writer(f).writerows(summary_rows)

print(f"\nDONE. Upload for the heatmap to be built:")
print(f"  {summary_path}")
for gene in regions:
    print(f"  {out_dir}/{gene}_signal_profile.csv")
