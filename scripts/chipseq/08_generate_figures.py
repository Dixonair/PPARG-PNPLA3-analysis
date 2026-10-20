#!/usr/bin/env python3
"""
Generate the two main ChIP-seq figures and the multi-gene heatmap for the
PNPLA3/PPARG ChIP-seq analysis.

INPUT (edit paths):
  - PNPLA3_signal_profile.csv, FABP4_signal_profile.csv   (per-bin signal, 10 samples)
  - PNPLA3_quantitative_signal_metrics.csv, FABP4_quantitative_signal_metrics.csv
  - bound_genes_promoter.csv
  - combined_8gene_heatmap_data.csv  (PNPLA3/FABP4 + ADIPOQ/CD36/PLIN1 + MYOD1/RHO/INS)
  - ASCRosi_DE_with_ChIP_annotation.csv

OUTPUT:
  - MAIN_Figure_A_tracks.pdf                  (full 10-sample signal tracks, PNPLA3 + FABP4)
  - MAIN_Figure_B_specificity_function.pdf    (quantitative/target-panel/heatmap/integration)

Coordinates used (hg38, verified against UCSC/GEO metadata during analysis):
  PNPLA3 TSS = chr22:43,923,805   FABP4 TSS (promoter peak centre) = chr8:81,483,350
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import fisher_exact

# Colour scheme used consistently across this analysis
COL_HASC = "#5B4B8A"
COL_TISSUE = "#E4A62B"
COL_BOUND = "#B565A7"
COL_UNBOUND = "#8AA29E"

SAMPLE_ORDER = ['GSM3177950', 'GSM3177951', 'GSM3177952', 'GSM3177953', 'GSM3177954',
                'GSM3416379', 'GSM3416380', 'GSM3416381', 'GSM3416382', 'GSM3416383']
PNPLA3_TSS = 43923805
FABP4_TSS = 81483350


def build_main_figure_A():
    fig = plt.figure(figsize=(16, 16.8))
    gs = gridspec.GridSpec(10, 4, hspace=0.15, wspace=0.4, top=0.96)

    pnpla3 = pd.read_csv('PNPLA3_signal_profile.csv')
    fabp4 = pd.read_csv('FABP4_signal_profile.csv')
    ymax_p = pnpla3[SAMPLE_ORDER].values.max()
    ymax_f = fabp4[SAMPLE_ORDER].values.max()

    for i, sample in enumerate(SAMPLE_ORDER):
        is_hasc = sample.startswith('GSM3177')
        color = COL_HASC if is_hasc else COL_TISSUE

        ax1 = fig.add_subplot(gs[i, 0:2])
        ax1.fill_between(pnpla3['position'], pnpla3[sample], color=color, alpha=0.85, linewidth=0)
        ax1.axvline(PNPLA3_TSS, color='firebrick', linestyle='--', linewidth=0.7, alpha=0.7)
        ax1.set_ylim(0, ymax_p * 1.05)
        ax1.set_yticks([])
        ax1.set_ylabel(sample.replace('GSM', ''), fontsize=7, rotation=0, ha='right', va='center')
        for sp in ['top', 'right', 'left']:
            ax1.spines[sp].set_visible(False)
        if i == 0:
            ax1.set_title('A. PPARG signal at PNPLA3 \u2014 every individual sample',
                           fontsize=11, loc='left', fontweight='bold')
        if i < 9:
            ax1.set_xticks([])
        else:
            ax1.set_xlabel('chr22 position (hg38)', fontsize=8)

        ax2 = fig.add_subplot(gs[i, 2:4])
        ax2.fill_between(fabp4['position'], fabp4[sample], color=color, alpha=0.85, linewidth=0)
        ax2.axvline(FABP4_TSS, color='firebrick', linestyle='--', linewidth=0.7, alpha=0.7)
        ax2.set_ylim(0, ymax_f * 1.05)
        ax2.set_yticks([])
        for sp in ['top', 'right', 'left']:
            ax2.spines[sp].set_visible(False)
        if i == 0:
            ax2.set_title('B. PPARG signal at FABP4 (positive control) \u2014 every individual sample',
                           fontsize=11, loc='left', fontweight='bold')
        if i < 9:
            ax2.set_xticks([])
        else:
            ax2.set_xlabel('chr8 position (hg38)', fontsize=8)

    fig.suptitle('Figure | PPARG ChIP-seq occupancy at PNPLA3 and FABP4 is reproducible '
                 'across all individual samples\n(purple = hASC, n=5; gold = adipose tissue, '
                 'n=5; dashed line = TSS; RPGC-normalised signal)', fontsize=12, y=0.985)
    plt.savefig('MAIN_Figure_A_tracks.pdf', bbox_inches='tight')
    plt.close()
    print("Saved MAIN_Figure_A_tracks.pdf")


def build_main_figure_B():
    fig = plt.figure(figsize=(15, 14))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1.1, 1], hspace=0.55, wspace=0.3)

    # Panel A: quantitative signal, PNPLA3 vs FABP4
    ax = fig.add_subplot(gs[0, 0])
    pnpla3_m = pd.read_csv('PNPLA3_quantitative_signal_metrics.csv')
    fabp4_m = pd.read_csv('FABP4_quantitative_signal_metrics.csv')
    for gi, (gene, df) in enumerate([('FABP4', fabp4_m), ('PNPLA3', pnpla3_m)]):
        for group, color, offset in [('hASC', COL_HASC, -0.15), ('tissue', COL_TISSUE, 0.15)]:
            vals = df[df.group == group]['mean_signal']
            x = np.full(len(vals), gi + offset) + np.random.normal(0, 0.02, len(vals))
            ax.scatter(x, vals, color=color, s=70, alpha=0.85, edgecolor='black', zorder=3,
                       label=group if gi == 0 else None)
            ax.errorbar([gi + offset], [vals.mean()], yerr=[vals.sem()], fmt='_', color='black',
                        markersize=16, capsize=4, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['FABP4', 'PNPLA3'], fontsize=10)
    ax.set_ylabel('Mean normalised signal\n(TSS\u00B12kb, RPGC)', fontsize=9)
    ax.set_title('A. Quantitative occupancy, every sample', fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=8)

    # Panel B: binding reproducibility vs canonical targets (PPARG shown separately)
    ax = fig.add_subplot(gs[0, 1])
    prom = pd.read_csv('bound_genes_promoter.csv')
    canonical_7 = ['ADIPOQ', 'CD36', 'FABP4', 'PCK1', 'PLIN1', 'LPL', 'SLC2A4']
    all_genes_shown = canonical_7 + ['PPARG', 'PNPLA3']
    sub = prom[prom['SYMBOL'].isin(all_genes_shown)]
    pivot_n = sub.pivot(index='SYMBOL', columns='dataset', values='n_samples_bound').reindex(all_genes_shown)
    x = np.arange(len(all_genes_shown)); width = 0.35
    ax.bar(x - width / 2, pivot_n['GSE115421_hASC'], width, label='hASC', color=COL_HASC)
    ax.bar(x + width / 2, pivot_n['GSE115421_tissue'], width, label='Tissue', color=COL_TISSUE)
    ax.axvline(6.5, color='black', linewidth=0.8, alpha=0.4)
    ax.axhline(3, color='grey', linestyle=':', linewidth=1)
    ax.set_xticks(x); ax.set_xticklabels(all_genes_shown, fontsize=8, rotation=30, ha='right')
    ax.get_xticklabels()[7].set_style('italic')
    ax.get_xticklabels()[8].set_fontweight('bold'); ax.get_xticklabels()[8].set_color('firebrick')
    ax.set_ylabel('Samples with peak (of 5)', fontsize=9)
    ax.set_title('B. Binding reproducibility vs canonical targets\n(PPARG shown separately)',
                 fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=8)

    # Panel C: multi-gene heatmap (positive + negative controls)
    ax = fig.add_subplot(gs[1, :])
    combined = pd.read_csv('combined_8gene_heatmap_data.csv')
    gene_order = ['ADIPOQ', 'CD36', 'FABP4', 'PLIN1', 'PNPLA3', 'MYOD1', 'RHO', 'INS']
    matrix = combined.pivot(index='gene', columns='sample', values='mean_signal').reindex(gene_order)[SAMPLE_ORDER]
    im = ax.imshow(np.log2(matrix + 1), cmap='viridis', aspect='auto')
    ax.set_xticks(range(10)); ax.set_xticklabels(SAMPLE_ORDER, rotation=90, fontsize=7)
    ax.set_yticks(range(8)); ax.set_yticklabels(gene_order, fontsize=10)
    ax.get_yticklabels()[4].set_color('firebrick'); ax.get_yticklabels()[4].set_fontweight('bold')
    ax.axhline(3.5, color='white', linewidth=1.5); ax.axhline(4.5, color='white', linewidth=1.5)
    ax.set_title('C. PNPLA3 signal vs established PPARG targets vs true negative-control loci '
                 '(MYOD1, RHO, INS)', fontsize=10, loc='left', fontweight='bold')
    plt.colorbar(im, ax=ax, fraction=0.025, label='log2(mean signal+1)')

    # Panel D/E: RNA-seq integration
    ax = fig.add_subplot(gs[2, 0])
    de = pd.read_csv('ASCRosi_DE_with_ChIP_annotation.csv')
    bound_fc = de[de.bound]['log2FC_Rosi_vs_DMSO']
    unbound_fc = de[~de.bound]['log2FC_Rosi_vs_DMSO']
    parts = ax.violinplot([unbound_fc, bound_fc], positions=[0, 1], showmeans=True)
    for pc, color in zip(parts['bodies'], [COL_UNBOUND, COL_BOUND]):
        pc.set_facecolor(color); pc.set_alpha(0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Not bound\n(n={len(unbound_fc)})', f'PPARG-bound\n(n={len(bound_fc)})'], fontsize=8)
    ax.set_ylabel('log2FC (Rosi vs DMSO)', fontsize=8)
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    pnpla3_fc = de[de.Gene == 'PNPLA3']['log2FC_Rosi_vs_DMSO'].values[0]
    ax.scatter([1], [pnpla3_fc], color='gold', s=110, edgecolor='black', zorder=5, marker='*')
    ax.annotate('PNPLA3', (1, pnpla3_fc), textcoords="offset points", xytext=(10, 0),
                fontweight='bold', fontsize=8)
    ax.set_title('D. Rosiglitazone response by\nbinding status', fontsize=10, loc='left', fontweight='bold')

    ax2 = fig.add_subplot(gs[2, 1])
    a = ((de.bound) & (de.induced)).sum(); b = ((de.bound) & (~de.induced)).sum()
    c = ((~de.bound) & (de.induced)).sum(); d = ((~de.bound) & (~de.induced)).sum()
    induced_pct = [100 * a / (a + b), 100 * c / (c + d)]
    ax2.bar(['PPARG-bound', 'Not bound'], induced_pct, color=[COL_BOUND, COL_UNBOUND], alpha=0.85)
    for i, v in enumerate(induced_pct):
        ax2.text(i, v + 0.08, f'{v:.1f}%', ha='center', fontweight='bold', fontsize=9)
    ax2.set_ylabel('% induced by rosiglitazone\n(padj<0.05)', fontsize=8)
    odds_ratio, p_fisher = fisher_exact([[a, b], [c, d]], alternative='greater')
    p_exp = int(np.floor(np.log10(p_fisher)))
    p_mant = p_fisher / (10 ** p_exp)
    ax2.set_title(f'E. Association: OR = {odds_ratio:.2f}\nFisher P = {p_mant:.1f}\u00d710$^{{{p_exp}}}$',
                  fontsize=10, loc='left', fontweight='bold')

    fig.suptitle('Figure | PNPLA3 occupancy is specific relative to true negative controls '
                 'and associated with rosiglitazone responsiveness', fontsize=12, y=0.995)
    plt.savefig('MAIN_Figure_B_specificity_function.pdf', bbox_inches='tight')
    plt.close()
    print("Saved MAIN_Figure_B_specificity_function.pdf")


if __name__ == "__main__":
    build_main_figure_A()
    build_main_figure_B()
    print("\nDone.")
