#!/usr/bin/env python3
"""
Extract raw bigWig signal values at PNPLA3 and FABP4 for locus-specific plotting.
Run from the project root (same directory as data/bigwig/).
"""
import pyBigWig
import csv
import os

bw_dir = os.path.join(os.getcwd(), "data/bigwig")
out_dir = os.path.join(os.getcwd(), "results/signal_extraction")
os.makedirs(out_dir, exist_ok=True)

samples = ["GSM3177950", "GSM3177951", "GSM3177952", "GSM3177953", "GSM3177954",
           "GSM3416379", "GSM3416380", "GSM3416381", "GSM3416382", "GSM3416383"]

# Windows based on the real peak coordinates already established in this project
regions = {
    "PNPLA3": ("chr22", 43921800, 43931400),   # TSS ~43,923,805, covers both peak clusters
    "FABP4":  ("chr8", 81482000, 81498200),    # covers all 4 known regulatory sites
}

for gene, region in regions.items():
    chrom, start, end = region
    bin_size = 20
    n_bins = (end - start) // bin_size
    positions = [start + i * bin_size for i in range(n_bins)]
    all_values = {}

    for s in samples:
        bw_path = os.path.join(bw_dir, f"{s}.bw")
        if not os.path.exists(bw_path):
            print(f"MISSING bigWig for {s}, filling with NA for {gene}")
            all_values[s] = [None] * n_bins
            continue
        bw = pyBigWig.open(bw_path)
        vals = bw.stats(chrom, start, end, nBins=n_bins, type="mean")
        all_values[s] = [v if v is not None else 0 for v in vals]
        bw.close()

    rows = [["position"] + samples]
    for i, pos in enumerate(positions):
        rows.append([pos] + [all_values[s][i] for s in samples])

    out_csv = os.path.join(out_dir, f"{gene}_signal_profile.csv")
    with open(out_csv, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"{gene}: wrote {out_csv} ({n_bins} bins x {len(samples)} samples)")

print("\nDone.")
