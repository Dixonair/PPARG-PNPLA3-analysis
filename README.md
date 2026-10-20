# PPARγ–PNPLA3 analysis code

Custom analysis code accompanying the manuscript *PPARγ engages proximal PPREs to regulate PNPLA3 in adipocytes and is associated with PNPLA3 I148M-mediated adipose tissue inflammation* ([journal / DOI to be added]).

Licence: MIT (see `LICENSE`). Archived version: [Zenodo DOI to be added].

## What is here

| Folder | Content |
|---|---|
| `scripts/chipseq/` | PPARG ChIP-seq reprocessing from raw reads (GSE115421), peak annotation, signal quantification, ChIP–RNA-seq integration, Hallmark enrichment, motif scan, ChIP figures (scripts 01–08), and the consolidated uniform pipeline `gse115421_pparg_chip_raw_reads.sh` |
| `scripts/rnaseq/` | RNA-seq/microarray analyses: hASC rosiglitazone (own data and GSE115421 RNA-seq), mouse adipocyte Pparg knockout (GSE62937), human subcutaneous WAT by PNPLA3 genotype, liver Pparg knockout/TZD (GSE162276) |
| `scripts/figures_tables/` | Assembly of the main Fig. 6 and Supplementary Fig. 6, Hallmark GSEA for the knockout and hASC datasets, the supplementary Hallmark tables (Word), and `make_source_data_GSE115421.py` (Source Data workbook for the ChIP-seq panels) |
| `scripts/chipseq_hepg2/` | PPARG ChIP-seq re-analysis of the public HepG2 data (ENCODE ENCSR130VQL): raw reads to peaks and tracks, window statistics at PNPLA3 and control genes, genome-wide percentiles, motif scan, integration with the liver and hASC RNA-seq tables, figure panels and Source Data (scripts 00–08, `hepg2_ppar_chip_raw_reads.sh`, `run_hepg2_all.sh`) |

Run order within `scripts/chipseq/` is numerical (01 → 08). Each script's header lists its inputs and outputs. The `rnaseq/` scripts are independent of each other.

## Uniform pipeline for all ten GSE115421 samples, and Source Data

`scripts/chipseq/gse115421_pparg_chip_raw_reads.sh` runs all ten samples (five hASC-derived adipocytes, GSM3177950–54; five adipose tissue, GSM3416379–83) through one pipeline with identical parameters, from SRA download to peaks, bigWig tracks, a QC table (`logs/qc_summary.tsv`) and a per-group >= 3-of-5 peak consensus. It uses the same commands as scripts 01 and 03 and adds a tool-version log. Test one sample first: `bash scripts/chipseq/gse115421_pparg_chip_raw_reads.sh 1`. Its header lists its limitations (no input control; 50 nt vs 100 nt reads; hg38 vs the authors' hg19 tracks).

`scripts/figures_tables/make_source_data_GSE115421.py [UP] [OUT] [xlsx]` writes the Source Data workbook (one sheet per panel, plus README, Files and Checks sheets) for the ChIP-seq panels (main Fig. 6e–l, Supplementary Figs. 7 and 8). `UP` is the folder with the signal profiles and peak tables, `OUT` the folder with the metrics, heatmap, RNA-seq and enrichment tables. Its README sheet lists the caveats (FABP4 window, z-score convention, paired t-test).

## HepG2 PPARG ChIP-seq (ENCODE ENCSR130VQL)

All scripts are in `scripts/chipseq_hepg2/` and use their own environment (`environment_hepg2.yml`; MACS3, in contrast to the MACS2 used for GSE115421). Raw reads to peaks, tracks and QC (two ChIP replicates, one input library; test with a small subset first):

```bash
conda env create -f environment_hepg2.yml && conda activate pparg_hepg2
INDEX=/path/to/bowtie2_index_GRCh38 CHIP=rep1.fastq.gz,rep2.fastq.gz CTRL=input.fastq.gz \
  bash scripts/chipseq_hepg2/hepg2_ppar_chip_raw_reads.sh
bash scripts/chipseq_hepg2/run_hepg2_all.sh      # downstream analysis, figure panels
OUT=results/hepg2_ppar GTF=ref/gencode.v46.basic.annotation.gtf.gz python3 scripts/chipseq_hepg2/08_export_source_data.py
```

`run_hepg2_all.sh` checks every input and tool first and stops with a message if something is missing. Inputs expected in `inputs/`, `ref/` or `data/`: the original GSE162276 and hASC RNA-seq tables (step 00 rebuilds the derived tables from them), GENCODE v46 annotation, the genome FASTA, and a Hallmark gene-set file (`hallmark_sets.pkl`, a pickled dict term -> human gene symbols). The header of each script lists its inputs and outputs. Steps 01 and 04 are superseded by `hepg2_ppar_chip_raw_reads.sh` and kept for the record.

## Data

No raw data are included. Public data used:
- GSE115421 (Hu et al. 2019): PPARG ChIP-seq (10 samples) and hASC rosiglitazone RNA-seq.
- ENCODE ENCSR130VQL (PPARG ChIP-seq, HepG2) with the input control ENCSR239QGH.
- GSE62937: adipocyte-specific inducible Pparg knockout microarray (GPL6887).
- GSE162276: hepatocyte-specific Pparg knockout and TZD RNA-seq.

Data generated for this study (hASC rosiglitazone/ligand experiments, human pioglitazone cohort, human adipose RNA-seq by PNPLA3 genotype, STAM mouse experiments): see the manuscript Data availability statement. Input file names expected by each script are given in its header; put them in `data/` (or set `PPARG_DATA`, `PPARG_RESULTS` for the R scripts, and `PPARG_DIR` for the Python scripts in `figures_tables/`).

## Software

See `environment.yml` (GSE115421 and the other analyses) and `environment_hepg2.yml` (HepG2 ChIP-seq) (Python/command-line tools) and the R package list in its comments. Record the exact versions used before release (`sessionInfo()`, `conda env export`).

## Quick start (ChIP-seq part)

```bash
conda env create -f environment.yml && conda activate pparg_pnpla3
# download Bowtie2 GRCh38 no-alt index into data/ref/GRCh38_noalt_as/ (see header of 01)
bash scripts/chipseq/01_reprocess_raw_reads.sh 1     # test with one sample, then without the argument
Rscript scripts/chipseq/02_annotate_peaks.R
bash scripts/chipseq/03_bam_to_bigwig.sh
python scripts/chipseq/04a_extract_locus_signal.py
python scripts/chipseq/04b_extract_multigene_signal.py
python scripts/chipseq/04c_compute_signal_metrics.py
python scripts/chipseq/05_integrate_with_rnaseq.py
python scripts/chipseq/06_pathway_enrichment.py
python scripts/chipseq/07_motif_scan.py
python scripts/chipseq/08_generate_figures.py
```

Figure/table assembly (all inputs and outputs in one folder):

```bash
export PPARG_DIR=/path/to/folder_with_inputs
python scripts/figures_tables/12_hallmark_gsea_mouse_and_hasc.py   # Hallmark GSEA (needs internet for gene sets once)
python scripts/figures_tables/09_rearranged_figure6.py             # needs figure_6.pdf (panels a–d artwork) in that folder
python scripts/figures_tables/10_make_hallmark_tables_docx.py       # supplementary Hallmark tables
python scripts/figures_tables/11_make_chip_hallmark_table_docx.py
python scripts/figures_tables/13_make_hallmark_tables_KO_hASC_docx.py
```
Scripts 10, 11 and 13 read their CSVs from the current directory; run them from `$PPARG_DIR`.

## Key results reproduced by this code

- PPARG peak at PNPLA3 (chr22:43,925,319–43,925,915, hg38; intron 1) called in 5/5 hASC and 4/5 adipose tissue samples; mean signal at the estimated upstream PPRE window (−790/−707) 0.9 versus 20.2 RPGC at the intron-1 peak.
- PPARG-bound genes are enriched among rosiglitazone-induced genes: 118/3,298 versus 85/13,372; OR 5.80, Fisher P = 1.3×10⁻³³.
- Hallmark adipogenesis is the most enriched set among PPARG-bound genes (FDR q = 6.9×10⁻²⁸).

## Limitations stated for transparency

1. ChIP peaks were called without an input control (none deposited); a 3-of-5 consensus is used for gene-level calls.
2. The ChIP–RNA-seq association (script 05) is correlative and uses a donor-averaged paired t-test, not a blocked model; do not quote it as the primary differential-expression result.
3. The motif scan (script 07) was inconclusive, including at its positive control; it is not used for any claim.
4. The signal window for FABP4 is truncated on the upstream side (TSS −1.35 kb to +2 kb) because extraction began at chr8:81,482,000.
5. Gene-set enrichment for the mouse knockout and hASC data (script 12) uses a variance-moderated t statistic (mouse) and the donor-paired t (hASC) in `gseapy` pre-ranked mode; results may differ slightly from `fgsea`/limma runs in the R scripts.
6. Pparg mRNA is not a valid knockout readout in GSE62937 and GSE162276 (a truncated transcript is detected; levels are higher, not lower).
7. `rnaseq/04_liver_Pparg_KO_and_TZD_GSE162276.R` contains `PATH/TO/...` placeholders for the count and metadata files; set them before running.
