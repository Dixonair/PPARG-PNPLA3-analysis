#!/usr/bin/env python3
"""
Integrate PPARG ChIP-seq binding calls (hASC) with matched rosiglitazone
RNA-sequencing data (GSE115421) to test whether PPARG-bound genes are
enriched among transcriptionally induced genes.

INPUT (edit paths below):
  - ASC_Rosi.csv              gene x sample expression matrix (GSE115421 RNA-seq)
                               sample columns named <Patient><Condition><Replicate>,
                               e.g. P1DMSO1, P1Rosi1 (5 patients, 3 reps/condition)
  - bound_genes_promoter.csv  output of 02_annotate_peaks.R

OUTPUT:
  - ASCRosi_DE_results_paired.csv       genome-wide paired DE results
  - ASCRosi_DE_with_ChIP_annotation.csv DE results annotated with ChIP binding status
  - integration_enrichment_stats.txt    Fisher's exact test result

Statistical approach: technical replicates are averaged to the patient level
before testing (5 patients, paired DMSO/Rosi), avoiding pseudoreplication from
treating 30 samples as 30 independent biological units. This is a Welch/paired-
t-test fallback; the manuscript's primary RNA-seq statistics should come from
limma with patient as a blocking factor (see Methods and the main hASC
reanalysis scripts), fit on all 30 samples directly.
"""
import pandas as pd
import numpy as np
import re
from scipy import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ---- paths (edit) -----------------------------------------------------------
ASC_ROSI_PATH = "ASC_Rosi.csv"
BOUND_GENES_PATH = "bound_genes_promoter.csv"
OUT_DIR = "."

# ---- 1. Load expression data, derive metadata from column names -------------
df = pd.read_csv(ASC_ROSI_PATH)
sample_cols = df.columns[1:].tolist()
gene_col = df.columns[0]

meta = []
for s in sample_cols:
    m = re.match(r'(P\d+)(DMSO|Rosi)(\d+)', s)
    if not m:
        raise ValueError(f"Sample name '{s}' does not match expected pattern "
                          f"<Patient><DMSO|Rosi><Replicate>, e.g. P1DMSO1")
    meta.append({'sample': s, 'patient': m.group(1), 'condition': m.group(2)})
meta = pd.DataFrame(meta)
print(f"Loaded {len(sample_cols)} samples: {meta['condition'].value_counts().to_dict()}, "
      f"{meta['patient'].nunique()} patients")

# ---- 2. Resolve duplicate gene symbols (keep highest-mean-expression row) ---
n_dup = df[gene_col].duplicated().sum()
print(f"Duplicated gene symbols: {n_dup} of {len(df)} rows")
df['_mean_expr'] = df[sample_cols].mean(axis=1)
df_dedup = df.sort_values('_mean_expr', ascending=False).drop_duplicates(subset=gene_col, keep='first')
df_indexed = df_dedup.set_index(gene_col)[sample_cols]

# ---- 3. Patient-level averaging (avoids pseudoreplication) + paired test ----
patients = sorted(meta['patient'].unique())
patient_means = {}
for p in patients:
    for cond in ['DMSO', 'Rosi']:
        cols = meta[(meta.patient == p) & (meta.condition == cond)]['sample'].tolist()
        patient_means[f"{p}_{cond}"] = df_indexed[cols].mean(axis=1)
patient_avg = pd.DataFrame(patient_means)
log_avg = np.log2(patient_avg + 1)

dmso_cols = [f"{p}_DMSO" for p in patients]
rosi_cols = [f"{p}_Rosi" for p in patients]
dmso_vals = log_avg[dmso_cols].values
rosi_vals = log_avg[rosi_cols].values

t_stats, pvals = stats.ttest_rel(rosi_vals, dmso_vals, axis=1)
log2fc = rosi_vals.mean(axis=1) - dmso_vals.mean(axis=1)

res_df = pd.DataFrame({'Gene': log_avg.index, 'log2FC_Rosi_vs_DMSO': log2fc,
                        't_stat': t_stats, 'pvalue': pvals}).dropna(subset=['pvalue'])
res_df['padj'] = multipletests(res_df['pvalue'], method='fdr_bh')[1]
res_df = res_df.sort_values('padj')
res_df.to_csv(f"{OUT_DIR}/ASCRosi_DE_results_paired.csv", index=False)
print(f"\nDE testing: {len(res_df)} genes, {(res_df['padj']<0.05).sum()} significant (padj<0.05)")

# ---- 4. Annotate with ChIP binding status and test enrichment ---------------
bound = pd.read_csv(BOUND_GENES_PATH)
bound_hasc = set(bound[bound['dataset'] == 'GSE115421_hASC']['SYMBOL'])
print(f"PPARG promoter-bound genes (hASC, consensus): {len(bound_hasc)}")

res_df['bound'] = res_df['Gene'].isin(bound_hasc)
res_df['induced'] = (res_df['padj'] < 0.05) & (res_df['log2FC_Rosi_vs_DMSO'] > 0)
res_df.to_csv(f"{OUT_DIR}/ASCRosi_DE_with_ChIP_annotation.csv", index=False)

a = ((res_df.bound) & (res_df.induced)).sum()
b = ((res_df.bound) & (~res_df.induced)).sum()
c = ((~res_df.bound) & (res_df.induced)).sum()
d = ((~res_df.bound) & (~res_df.induced)).sum()
odds_ratio, pval_fisher = fisher_exact([[a, b], [c, d]], alternative='greater')

with open(f"{OUT_DIR}/integration_enrichment_stats.txt", "w") as f:
    f.write("PPARG-bound vs rosiglitazone-induced gene enrichment\n")
    f.write("=" * 55 + "\n")
    f.write(f"                 Induced   Not induced\n")
    f.write(f"PPARG-bound      {a:6d}    {b:6d}\n")
    f.write(f"Not bound        {c:6d}    {d:6d}\n\n")
    f.write(f"Odds ratio: {odds_ratio:.2f}\n")
    f.write(f"Fisher's exact P (one-sided): {pval_fisher:.4g}\n")
    f.write(f"% bound genes induced: {100*a/(a+b):.1f}%\n")
    f.write(f"% unbound genes induced: {100*c/(c+d):.1f}%\n")

print(f"\nOdds ratio: {odds_ratio:.2f}, Fisher P = {pval_fisher:.4g}")
pnpla3_row = res_df[res_df.Gene == 'PNPLA3']
if len(pnpla3_row):
    print(f"\nPNPLA3: log2FC={pnpla3_row.iloc[0]['log2FC_Rosi_vs_DMSO']:.3f}, "
          f"padj={pnpla3_row.iloc[0]['padj']:.3f}, bound={pnpla3_row.iloc[0]['bound']}")
print("\nDone. Outputs written to", OUT_DIR)
