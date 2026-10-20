# --- Portable paths: set these two lines (or run from the project root with data/ and results/ folders) ---
DATA_DIR <- Sys.getenv("PPARG_DATA", "data"); RESULTS_DIR <- Sys.getenv("PPARG_RESULTS", "results")

#!/usr/bin/env Rscript
# ==============================================================================
# GSE162276: hepatic RNA-seq, Pparg-Hep-KO vs WT, +/- rosiglitazone (Lee et al. 2021)
# Two comparisons, matching your existing figure:
#   A) HFCF ROSI vs HFCF CON   (rosiglitazone effect, Pparg-intact liver)
#   B) HFCF-KO   vs HFCF-WT    (genotype effect)
#
# *** UNTESTED against real data -- GEO is unreachable from the sandbox that
# *** wrote this script (same restriction hit throughout the ChIP-seq work).
# *** Syntax-checked only. Run Section 2's diagnostics FIRST and read them
# *** before trusting anything downstream -- do not assume column names.
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(msigdbr); library(clusterProfiler)
  library(org.Mm.eg.db); library(pathview); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(pheatmap)
})

## ---- 1. Paths (EDIT) ---------------------------------------------------------
count_path <- "PATH/TO/your_countmatrix.csv"
meta_path  <- "PATH/TO/your_metadata.csv"
out_dir    <- "results/GSE162276"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 2. Load + DIAGNOSE (read this output before proceeding) ----------------
counts_raw <- read.csv(count_path, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)
cat("Count matrix:", nrow(counts_raw), "x", ncol(counts_raw), "\n")
cat("First col (should be gene symbol):", colnames(counts_raw)[1], "\n")
print(head(counts_raw[[1]]))
cat("\nMetadata:\n"); print(meta)
cat("\nUnique values in metadata condition column(s):\n")
for (col in colnames(meta)) if (is.character(meta[[col]])) { cat(col, ":\n"); print(unique(meta[[col]])) }

# Raw-vs-normalised diagnostic (this exact ambiguity caused real problems earlier
# in this project for other datasets -- check it explicitly, don't assume).
num_check <- as.numeric(unlist(counts_raw[1:min(5, nrow(counts_raw)), sample_cols_check <- setdiff(colnames(counts_raw), colnames(counts_raw)[1])]))
cat("\nSample values (first 5 genes, all samples):\n"); print(head(num_check, 15))
cat("All integer-like?", all(num_check == round(num_check), na.rm = TRUE), "\n")
libsize <- colSums(counts_raw[, sample_cols_check])
cat("Library size CV (variable=likely raw counts; near-zero=likely pre-normalised):",
    round(sd(libsize)/mean(libsize), 4), "\n")
cat("*** If values are NOT integer-like and/or library-size CV is near zero, this ***\n")
cat("*** is pre-normalised data (as GSE115421 and the human WAT dataset were) --  ***\n")
cat("*** DESeq2 in Section 4 is then WRONG; switch to limma on log2 values,       ***\n")
cat("*** matching the pattern already used for GSE115421 in this project.        ***\n\n")

## ---- 3. EDIT to match your actual metadata -----------------------------------
SAMPLE_COL <- "sample"       # <-- EDIT
GROUP_COL  <- "condition"    # <-- EDIT
# Exact strings as they appear in your metadata (from Section 2 output):
ROSI_LEVEL    <- "HFCF ROSI"
ROSI_CONTROL  <- "HFCF CON"
KO_LEVEL      <- "HFCF-KO"
WT_LEVEL      <- "HFCF-WT"

colnames(meta)[colnames(meta) == SAMPLE_COL] <- "sample"
colnames(meta)[colnames(meta) == GROUP_COL]  <- "condition"

gene_col <- colnames(counts_raw)[1]
n_dup <- sum(duplicated(counts_raw[[gene_col]]))
cat("\nDuplicate gene symbols:", n_dup, "\n")
if (n_dup > 0) {
  cat("Aggregating (SUM) duplicate symbols -- valid for RNA-seq counts.\n")
  counts_raw <- counts_raw %>% group_by(.data[[gene_col]]) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
}

sample_cols <- setdiff(colnames(counts_raw), gene_col)
counts_all <- counts_raw %>% column_to_rownames(gene_col)

## ---- 4. Run BOTH comparisons via a reusable function -------------------------

run_comparison <- function(counts_all, meta, group_col, level_a, level_b, label) {
  cat("\n\n========== COMPARISON:", label, "(", level_b, "vs", level_a, ") ==========\n")
  sub_meta <- meta %>% filter(.data[[group_col]] %in% c(level_a, level_b))
  sub_meta[[group_col]] <- factor(sub_meta[[group_col]], levels = c(level_a, level_b))
  common <- intersect(sub_meta$sample, colnames(counts_all))
  cat(length(common), "of", nrow(sub_meta), "metadata samples matched to count columns\n")
  sub_meta <- sub_meta[match(common, sub_meta$sample), ]
  cts <- round(as.matrix(counts_all[, common]))
  mode(cts) <- "integer"
  cts[cts < 0] <- 0  # guard: if the file is normalised counts, not raw, this will be wrong -- check Section 2 diagnostics

  keep <- rowSums(cts >= 10) >= min(table(sub_meta[[group_col]]))
  cts <- cts[keep, ]
  cat("Genes retained after filtering:", nrow(cts), "\n")

  dds <- DESeqDataSetFromMatrix(countData = cts, colData = sub_meta,
                                 design = as.formula(paste0("~", group_col)))
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(group_col, level_b, level_a)) %>%
    as.data.frame() %>% rownames_to_column("Gene") %>% arrange(padj)
  write.csv(res, file.path(out_dir, paste0("DESeq2_", label, ".csv")), row.names = FALSE)
  cat("N significant (FDR<0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")

  # GSEA (Hallmark + GO), ranked by Wald statistic
  ranked <- res %>% filter(!is.na(stat)) %>% group_by(Gene) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(desc(stat))
  rank_vec <- setNames(ranked$stat, ranked$Gene)

  hallmark <- msigdbr(species = "Mus musculus", collection = "MH", db_species = "MM") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  gsea_h <- fgsea(pathways = hallmark, stats = rank_vec, minSize = 10, maxSize = 500) %>% arrange(padj)
  write.csv(as.data.frame(gsea_h %>% select(-leadingEdge)),
            file.path(out_dir, paste0("GSEA_Hallmark_", label, ".csv")), row.names = FALSE)
  cat("Top Hallmark pathways:\n"); print(head(gsea_h[, c("pathway","NES","padj")], 10))

  # PCA
  vsd <- vst(dds, blind = FALSE)
  pca <- prcomp(t(assay(vsd)[order(-apply(assay(vsd),1,var))[1:min(2000,nrow(vsd))],]), scale.=TRUE)
  pca_df <- as.data.frame(pca$x[,1:2]) %>% rownames_to_column("sample") %>% left_join(sub_meta, by="sample")
  p <- ggplot(pca_df, aes(PC1, PC2, color = .data[[group_col]])) + geom_point(size=4) +
    theme_bw() + labs(title = paste("PCA:", label))
  ggsave(file.path(out_dir, paste0("PCA_", label, ".pdf")), p, width=6, height=5)

  # Curated panel, matching your existing figure
  panel_genes <- c("Cd36","Fabp4","Pnpla2","Pnpla3","Pparg","Slc2a4","Cebpa","Fasn","Nr1h3","Ppargc1a","Scd1","Srebf1")
  panel_present <- intersect(panel_genes, rownames(assay(vsd)))
  cat("Panel genes present:", length(panel_present), "of", length(panel_genes), "\n")
  hm <- assay(vsd)[panel_present, common]
  hm_z <- t(scale(t(hm)))
  pdf(file.path(out_dir, paste0("panel_heatmap_", label, ".pdf")), width=7, height=6)
  pheatmap(hm_z, cluster_cols = FALSE, main = paste("Panel:", label))
  dev.off()

  list(dds = dds, res = res, gsea = gsea_h)
}

result_rosi <- run_comparison(counts_all, meta, "condition", ROSI_CONTROL, ROSI_LEVEL, "Rosi_vs_Control")
result_ko   <- run_comparison(counts_all, meta, "condition", WT_LEVEL, KO_LEVEL, "KO_vs_WT")

## ---- 5. GO enrichment (ORA) on nominal-p genes, KO vs WT (mirrors human WAT) -
sig_genes_ko <- result_ko$res %>% filter(pvalue < 0.01) %>% pull(Gene)
ego <- enrichGO(gene = sig_genes_ko, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ego), file.path(out_dir, "GO_ORA_KO_vs_WT.csv"), row.names = FALSE)

cat("\n\n===== DONE =====\nOutputs in:", out_dir, "\n")
cat("REMINDER: this ran on Section 2's untested diagnostics -- if column names,\n")
cat("group level strings, or the raw-vs-normalised-counts assumption didn't match\n")
cat("your actual file, results above are wrong. Check the printed diagnostics.\n")
