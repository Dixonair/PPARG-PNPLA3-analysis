#!/usr/bin/env Rscript
# ==============================================================================
# STEP 1 (run on YOUR machine): peak files -> bound-gene tables + PNPLA3 check
# Datasets: GSE115421 (Hu et al. 2019; hASC-derived adipocytes, 5 patients;
#           PPARG ChIP-seq, hg38)  and  GSE41578 (Galhardo et al.; SGBS
#           adipocytes; PPARG, CEBPA, LXR, H3K4me3; hg19)
#
# Requires Bioconductor (see header of the package README).
#
# WHAT YOU MUST SUPPLY: one peak file per sample (BED or narrowPeak; first three
# columns chr/start/end). Check the GEO "Supplementary files" for each series.
# GSE115421 is raw reads only in the SRA table, and the table lists NO input/IgG
# runs (ten PPARG ChIP runs, nothing else). Check GEO for input samples; if there
# are none, peaks must be called without a control (MACS2 with no -c), which is
# less specific -- state this in Methods and lean on the multi-sample consensus.
# GSE41578 does include an input sample (Galhardo et al. data descriptor).
#
# If you re-align GSE41578 to hg38 yourself, set its `build` entries below to "hg38".
# Genome builds MUST match the TxDb used. GSE115421 = hg38 (per your earlier
# work), GSE41578 = hg19. Confirm both on the GEO pages before trusting output.
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(ChIPseeker); library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(dplyr)
})

out_dir <- "results/step1_peaks_to_genes"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- EDIT: one row per peak file --------------------------------------------
peaks <- data.frame(
  # GSE115421 has TEN PPARG ChIP samples (SraRunTable, one BioProject PRJNA474921):
  #   GSM3177950-54 = hASC-derived adipocytes, 5 patients, 50 bp SE, agent "none"
  #   GSM3416379-83 = human ADIPOSE TISSUE, 5 samples, 100 bp SE (donor/depot
  #                   details are NOT in the run table -- read them off GEO)
  # GSM -> patient (P1..P5) mapping is also not in the table; GSM ids are used.
  # `pattern` is matched (case-insensitive) against file names inside peak_dir,
  # so MACS2 output like GSM3177950_peaks.narrowPeak is found without renaming.
  pattern = c(paste0("GSM", 3177950:3177954), paste0("GSM", 3416379:3416383),
              "GSE41578.*PPARG", "GSE41578.*CEBPA", "GSE41578.*LXR", "GSE41578.*H3K4me3"),
  dataset = c(rep("GSE115421_hASC", 5), rep("GSE115421_tissue", 5), rep("GSE41578", 4)),
  sample  = c(paste0("GSM", 3177950:3177954), paste0("GSM", 3416379:3416383),
              "SGBS","SGBS","SGBS","SGBS"),
  factor  = c(rep("PPARG", 10), "PPARG","CEBPA","LXR","H3K4me3"),
  build   = c(rep("hg38", 10), rep("hg19", 4)),
  stringsAsFactors = FALSE
)
peak_dir <- "data/peaks"      # folder holding the peak files (relative to getwd() = project root)
PROMOTER_WINDOW <- 3000      # bp either side of TSS = "promoter-bound"
NEAR_WINDOW     <- 50000     # bp = "nearby (promoter or enhancer)" assignment
MIN_DONORS      <- 3         # GSE115421 (each group of 5): bound if seen in >= 3 of 5

txdbs <- list(hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene,
              hg19 = TxDb.Hsapiens.UCSC.hg19.knownGene)

read_peaks <- function(f) {
  d <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  # Standard MACS2 narrowPeak has 10 columns: chrom, start, end, name, score,
  # strand, signalValue, pValue, qValue, summit (offset from start). The
  # original version of this function only read columns 1-3, silently
  # discarding signalValue/pValue/qValue -- needed for any analysis relating
  # continuous peak strength to another variable (e.g. RNA-seq response).
  gr <- GRanges(d[[1]], IRanges(d[[2]] + 1, d[[3]]))   # BED is 0-based
  if (ncol(d) >= 10) {
    mcols(gr)$signalValue <- d[[7]]
    mcols(gr)$pValue_peak <- d[[8]]   # suffixed _peak to avoid clashing with
    mcols(gr)$qValue_peak <- d[[9]]   # annotatePeak's own output columns
    mcols(gr)$summit <- d[[2]] + d[[10]]  # absolute summit position
  } else {
    warning(basename(f), " has only ", ncol(d), " columns -- not a standard ",
            "10-column narrowPeak file. signalValue/pValue/qValue/summit will ",
            "be NA for this sample.")
  }
  gr
}

cat("Working directory:", getwd(), "\n")
cat("Looking for peak files in:", normalizePath(peak_dir, mustWork = FALSE), "\n")
if (!dir.exists(peak_dir)) stop("Folder '", peak_dir, "' does not exist here. Create it and put the ",
                                "peak files in it, or set peak_dir / setwd() correctly.")
cat("Files found there:\n"); print(list.files(peak_dir))

annotated <- lapply(seq_len(nrow(peaks)), function(i) {
  p <- peaks[i, ]
  f <- list.files(peak_dir, pattern = p$pattern, full.names = TRUE, ignore.case = TRUE)
  f <- f[!grepl("_track", f, ignore.case = TRUE)]  # exclude UCSC-upload files (have a
                                                     # non-BED "track type=..." header line
                                                     # prepended; not raw peak data)
  np <- grep("narrowPeak$", f, value = TRUE)          # prefer MACS2 narrowPeak; ignore .xls / _summits.bed
  if (length(np) > 0) f <- np else f <- f[!grepl("summits|\\.xls$", f)]
  if (length(f) == 0) { message("MISSING (no file matching '", p$pattern, "')"); return(NULL) }
  if (length(f) > 1) message("Several files match '", p$pattern, "'; using ", basename(f[1]))
  gr <- read_peaks(f[1])
  ann <- annotatePeak(gr, TxDb = txdbs[[p$build]], tssRegion = c(-PROMOTER_WINDOW, PROMOTER_WINDOW),
                      annoDb = "org.Hs.eg.db", verbose = FALSE)
  df <- as.data.frame(ann)
  df$dataset <- p$dataset; df$sample <- p$sample; df$factor <- p$factor; df$build <- p$build
  df$file <- basename(f[1])
  df
})
n_found <- sum(!vapply(annotated, is.null, logical(1)))
cat(n_found, "of", nrow(peaks), "peak files found and annotated\n")

# Verify the new columns actually survived annotatePeak() -- this cannot be
# tested without Bioconductor access, so check it explicitly on first run
# rather than silently trust it.
first_real <- Filter(Negate(is.null), annotated)
if (length(first_real) > 0) {
  new_cols <- c("signalValue", "pValue_peak", "qValue_peak", "summit")
  present <- new_cols %in% colnames(first_real[[1]])
  cat("\nsignalValue/pValue/qValue/summit columns present after annotatePeak()?\n")
  print(setNames(present, new_cols))
  if (!all(present)) {
    cat("*** WARNING: some columns were dropped by annotatePeak(). If this  ***\n")
    cat("*** happens, extract signalValue from the GRanges BEFORE calling  ***\n")
    cat("*** annotatePeak() and join it back in afterwards by peak coord-  ***\n")
    cat("*** inate instead of relying on it surviving annotatePeak().      ***\n")
  }
}

if (n_found == 0) stop("No peak files were found, so there is nothing to annotate. This is a missing-input ",
                       "problem, not a code error: put peak files (BED/narrowPeak) in '", peak_dir, "'.")
if (n_found < nrow(peaks)) warning("Some datasets are missing; consensus counts below use only what was found.")
annotated <- dplyr::bind_rows(annotated)
write.csv(annotated, file.path(out_dir, "all_peaks_annotated.csv"), row.names = FALSE)

## ---- bound-gene tables --------------------------------------------------------
bound <- annotated %>% filter(!is.na(SYMBOL)) %>%
  mutate(abs_dist = abs(distanceToTSS),
         promoter = abs_dist <= PROMOTER_WINDOW, near = abs_dist <= NEAR_WINDOW)

summarise_bound <- function(df, window_col) {
  df %>% filter(.data[[window_col]]) %>%
    group_by(dataset, factor, SYMBOL) %>%
    summarise(n_samples_bound = n_distinct(sample), min_abs_dist = min(abs_dist), .groups = "drop")
}
for (w in c("promoter", "near")) {
  tab <- summarise_bound(bound, w)
  # GSE115421 consensus: keep genes seen in >= MIN_DONORS donors
  tab <- tab %>% filter(!(startsWith(dataset, "GSE115421") & n_samples_bound < MIN_DONORS))
  write.csv(tab, file.path(out_dir, paste0("bound_genes_", w, ".csv")), row.names = FALSE)
  cat(w, "-bound genes per dataset/factor:\n", sep = ""); print(table(tab$dataset, tab$factor))
}

## ---- PNPLA3-specific check ------------------------------------------------------
# Reports every peak within 10 kb of a PNPLA3 TSS in every sample, so you can
# confirm the peak sits at PNPLA3 and not a neighbouring gene (the earlier
# EFCAB6 mix-up is exactly what this guards against) and compare its position
# with your two PPREs (-713/-707 and -790/-780 from the TSS).
pnpla3 <- annotated %>% filter(SYMBOL == "PNPLA3", abs(distanceToTSS) <= 10000) %>%
  select(dataset, sample, factor, build, seqnames, start, end, distanceToTSS, annotation)
write.csv(pnpla3, file.path(out_dir, "PNPLA3_peaks.csv"), row.names = FALSE)
cat("\nPNPLA3 peaks within 10 kb of a TSS:\n"); print(pnpla3)
for (ds in c("GSE115421_hASC", "GSE115421_tissue")) {
  cat("\nSamples with a PNPLA3 promoter peak (", ds, "): ",
      length(unique(pnpla3$sample[pnpla3$dataset == ds & abs(pnpla3$distanceToTSS) <= PROMOTER_WINDOW])),
      " of 5\n", sep = "")
}
cat("\nDONE. Upload results/step1_peaks_to_genes/bound_genes_promoter.csv, bound_genes_near.csv and\n",
    "PNPLA3_peaks.csv (small files) so step 2 (integration with the rosiglitazone RNA-seq,\n",
    "GSEA, GO and network analysis) can be run against real data.\n", sep = "")
