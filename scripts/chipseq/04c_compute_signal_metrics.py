#!/usr/bin/env python3
"""
Compute per-sample quantitative PPARG occupancy metrics with ONE consistent
window definition (TSS +/- 2 kb) for every gene, and assemble the 8-gene
heatmap input table.

INPUT (outputs of 04a and 04b):
  - PNPLA3_signal_profile.csv, FABP4_signal_profile.csv   (04a; wider windows)
  - multigene_heatmap_summary.csv                          (04b; ADIPOQ, CD36,
        PLIN1, MYOD1, RHO, INS -- already extracted as TSS +/- 2 kb)

OUTPUT:
  - PNPLA3_quantitative_signal_metrics.csv
  - FABP4_quantitative_signal_metrics.csv
  - combined_8gene_heatmap_data.csv

WINDOW NOTES (reported, not hidden):
  * PNPLA3 TSS = chr22:43,923,805 (hg38). The TSS +/- 2 kb window contains
    ~91% of the summed signal of the intron-1 peak (chr22:43,925,319-43,925,915);
    the peak's right edge extends slightly beyond the window.
  * FABP4 marker = chr8:81,483,350 (approximate TSS, taken from the centre of
    the promoter-proximal peak). The 04a extraction began at chr8:81,482,000,
    so the FABP4 window is truncated on the upstream side (TSS-1.35 kb to
    TSS+2 kb, 168 of 200 bins). To obtain the full window, set the FABP4
    region in 04a to start at 81,481,300 and re-run.
  * This window covers only the FABP4 promoter-proximal peak, not the three
    downstream peaks shown in the tracks.
"""
import pandas as pd

SAMPLES = ["GSM3177950", "GSM3177951", "GSM3177952", "GSM3177953", "GSM3177954",
           "GSM3416379", "GSM3416380", "GSM3416381", "GSM3416382", "GSM3416383"]
TSS = {"PNPLA3": 43923805, "FABP4": 81483350}
HALF_WINDOW = 2000
BIN_SIZE = 20


def group_of(sample):
    return "hASC" if sample.startswith("GSM3177") else "tissue"


combined_rows = []
for gene, tss in TSS.items():
    prof = pd.read_csv(f"{gene}_signal_profile.csv")
    win = prof[(prof["position"] >= tss - HALF_WINDOW) & (prof["position"] <= tss + HALF_WINDOW)]
    n_expected = 2 * HALF_WINDOW // BIN_SIZE
    print(f"{gene}: {len(win)} bins in TSS+/-{HALF_WINDOW} bp window (full window = {n_expected})"
          + ("  ** TRUNCATED **" if len(win) < n_expected else ""))
    rows = []
    for s in SAMPLES:
        vals = win[s].values
        rows.append({"sample": s, "group": group_of(s), "mean_signal": vals.mean(),
                     "max_signal": vals.max(), "auc": vals.sum() * BIN_SIZE})
        combined_rows.append({"gene": gene, "sample": s, "group": group_of(s),
                              "mean_signal": vals.mean()})
    pd.DataFrame(rows).to_csv(f"{gene}_quantitative_signal_metrics.csv", index=False)

others = pd.read_csv("multigene_heatmap_summary.csv")
combined = pd.concat([pd.DataFrame(combined_rows), others], ignore_index=True)
combined.to_csv("combined_8gene_heatmap_data.csv", index=False)
print(f"\ncombined_8gene_heatmap_data.csv: {combined['gene'].nunique()} genes x "
      f"{combined['sample'].nunique()} samples = {len(combined)} rows")
print(combined.groupby("gene")["mean_signal"].mean().round(2).to_string())
