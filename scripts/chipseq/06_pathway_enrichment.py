#!/usr/bin/env python3
"""
Exploratory, genome-wide Hallmark pathway enrichment among PPARG-bound genes
(hypergeometric test), using expressed genes in the matched RNA-seq dataset
as the statistical background.

INPUT (edit paths):
  - bound_genes_promoter.csv       output of 02_annotate_peaks.R
  - ASCRosi_DE_with_ChIP_annotation.csv   output of 05_integrate_with_rnaseq.py
                                     (used only to define the expressed-gene background)
  - hallmark_sets.pkl               MSigDB Hallmark gene sets, human symbols
                                     (see note below on sourcing this file)

OUTPUT:
  - ChIP_bound_genes_Hallmark_ORA.csv

NOTE ON GENE SETS: this script expects a pickled Python dict {term_name: [gene
symbols]} for the 50 MSigDB Hallmark sets. Obtain this via msigdbr in R
(species="Homo sapiens", collection="H") and export to the same format, or
via a Python MSigDB client. Do not substitute an unverified source.
"""
import pandas as pd
import numpy as np
import pickle
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

BOUND_GENES_PATH = "bound_genes_promoter.csv"
DE_ANNOTATED_PATH = "ASCRosi_DE_with_ChIP_annotation.csv"
HALLMARK_PKL_PATH = "hallmark_sets.pkl"
OUT_PATH = "ChIP_bound_genes_Hallmark_ORA.csv"

bound = pd.read_csv(BOUND_GENES_PATH)
bound_hasc = set(bound[bound['dataset'] == 'GSE115421_hASC']['SYMBOL'])

de = pd.read_csv(DE_ANNOTATED_PATH)
background = set(de['Gene'])  # expressed genes, not whole genome
bound_in_bg = bound_hasc & background
print(f"Bound genes: {len(bound_hasc)}, present in expression background: {len(bound_in_bg)}")

with open(HALLMARK_PKL_PATH, 'rb') as f:
    hallmark = pickle.load(f)

N = len(background)
n_bound = len(bound_in_bg)
results = []
for term, genes in hallmark.items():
    term_in_bg = set(genes) & background
    K = len(term_in_bg)
    if K < 5:
        continue
    overlap = term_in_bg & bound_in_bg
    k = len(overlap)
    pval = hypergeom.sf(k - 1, N, K, n_bound)
    results.append({'Term': term, 'K_pathway_genes': K, 'k_bound_and_in_pathway': k,
                     'expected_if_random': K * n_bound / N, 'pval': pval,
                     'genes': ';'.join(sorted(overlap))})

res_df = pd.DataFrame(results)
res_df['padj'] = multipletests(res_df['pval'], method='fdr_bh')[1]
res_df = res_df.sort_values('padj')
res_df.to_csv(OUT_PATH, index=False)

print(f"\nTop 10 enriched Hallmark pathways:")
print(res_df[['Term', 'K_pathway_genes', 'k_bound_and_in_pathway', 'padj']].head(10).to_string())
print(f"\nDone. Output: {OUT_PATH}")
