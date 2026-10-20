# --- Portable paths: set these two lines (or run from the project root with data/ and results/ folders) ---
DATA_DIR <- Sys.getenv("PPARG_DATA", "data"); RESULTS_DIR <- Sys.getenv("PPARG_RESULTS", "results")

#!/usr/bin/env Rscript
# ==============================================================================
# PPARG agonism (rosiglitazone) -> PNPLA3 and the adipocyte transcriptional
# program, in human patient-derived adipose stem cell (ASC) adipocytes.
#
# Dataset: Hu, Jiang, Guan et al. 2019, Cell Stem Cell 24(2):299-308.e296.
#          "Patient Adipose Stem Cell-Derived Adipocytes Reveal Genetic
#          Variation that Predicts Antidiabetic Drug Response."
#
# *** THIS SCRIPT HAS NOT BEEN TESTED AGAINST YOUR ACTUAL DATA. ***
# Unlike the GSE62937 (Wang et al.) pipeline, these files were not uploaded
# for verification, so every assumption below (column names, design, data
# type) is a best guess and is checked diagnostically at runtime rather than
# hard-coded. READ THE CONSOLE OUTPUT of Sections 2-4 carefully before
# trusting anything downstream -- if the printed diagnostics don't match what
# you expect, STOP and fix the EDIT-marked variables before proceeding.
#
# KEY DESIGN DIFFERENCES FROM THE WANG ET AL. (MOUSE KO) PIPELINE:
#  - Human data: gene symbols are upper-case (PPARG, PNPLA3, ADIPOQ, ...).
#  - This is rosiglitazone TREATMENT of patient-derived cells, not genetic
#    knockout -- so the expected direction is OPPOSITE to the KO study:
#    PPARG agonism should INCREASE Pnpla3 and other PPARG targets, not
#    decrease them. Do not reuse the KO study's directionality assumptions.
#  - This study profiles many patient lines to study genetic variation in
#    drug response, so the design may be PAIRED (same patient, Rosi vs
#    Vehicle) with an n much larger than the 3-vs-3 mouse study. This
#    script auto-detects a subject/patient ID column and switches to a
#    paired analysis if one is found; check the printed message.
#  - The expression matrix may be raw RNA-seq counts (unlike the mouse study,
#    which was normalised microarray intensity) -- this script checks for
#    integer-like values and variable library sizes, and will tell you which
#    branch (DESeq2-style raw counts vs limma on normalised values) is
#    appropriate. Do not assume one or the other; read the diagnostic output.
#
# VISUAL STYLE: deliberately distinct from the Wang et al. (mouse KO) figures
# -- teal/coral palette instead of blue/orange, violin+jitter instead of
# boxplots, a different network layout aesthetic -- so the two datasets'
# figures are never confusable at a glance in the manuscript.
# ==============================================================================

## ---- 0. Setup ---------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggrepel", "pheatmap")
installed <- rownames(installed.packages())
missing_pkgs <- setdiff(required_pkgs, installed)
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap)
})

# Distinct visual identity for this dataset (vs Wang et al. blue/orange)
COL_VEHICLE <- "#1B9E77"   # teal
COL_ROSI    <- "#D95F02"   # coral/burnt orange
PALETTE_DIV <- c("#2C7FB8", "white", "#D95F02")  # blue-white-coral for heatmap

## ---- 1. Paths (EDIT IF NEEDED) ----------------------------------------------

expr_path    <- file.path(DATA_DIR, "ASC_Rosi.csv")
metadata_path <- file.path(DATA_DIR, "metadata.csv")
out_dir      <- file.path(RESULTS_DIR, "hASC_rosi")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 2. Load and DIAGNOSE the expression matrix -----------------------------
## Do not skip reading this output. It tells you what kind of data you have.

expr_raw <- read.csv(expr_path, check.names = FALSE)
cat("Expression file dimensions:", nrow(expr_raw), "rows x", ncol(expr_raw), "columns\n")
cat("First column name (should be a gene ID/symbol column):", colnames(expr_raw)[1], "\n")
cat("First 5 column names:\n"); print(head(colnames(expr_raw), 5))
cat("First 5 row values in column 1 (gene symbols/IDs):\n"); print(head(expr_raw[[1]], 5))

sample_cols <- setdiff(colnames(expr_raw), colnames(expr_raw)[1])
gene_col <- colnames(expr_raw)[1]

numeric_check <- suppressWarnings(as.numeric(unlist(expr_raw[1:min(5, nrow(expr_raw)), sample_cols])))
cat("\nSample of raw values (first 5 genes, all samples):\n"); print(summary(numeric_check))
cat("Are these integer-like (consistent with raw RNA-seq counts)?\n")
print(all(numeric_check == round(numeric_check), na.rm = TRUE))

col_sums <- colSums(sapply(expr_raw[sample_cols], as.numeric), na.rm = TRUE)
cat("\nPer-sample column sums (library size proxy) -- if these vary substantially\n")
cat("across samples, this is likely RAW RNA-seq counts (use DESeq2-style approach,\n")
cat("Section 5A). If they are all nearly identical, this is likely PRE-NORMALISED\n")
cat("data (use limma directly on log-values, Section 5B), as in the Wang et al. study.\n")
print(summary(col_sums))
cat("Coefficient of variation of library sizes:", round(sd(col_sums)/mean(col_sums), 4), "\n")
cat("(CV < ~0.02 suggests pre-normalised/equalised data, as seen in GSE62937;\n")
cat(" CV > ~0.1 suggests raw or minimally-normalised counts with real depth variation.)\n\n")

## ---- 3. Load and DIAGNOSE the metadata ---------------------------------------

meta <- read.csv(metadata_path, check.names = FALSE)
cat("Metadata dimensions:", nrow(meta), "rows x", ncol(meta), "columns\n")
cat("Metadata columns:\n"); print(colnames(meta))
cat("First rows:\n"); print(head(meta))
for (col in colnames(meta)) {
  if (is.character(meta[[col]]) || length(unique(meta[[col]])) < 20) {
    cat("\nUnique values in '", col, "':\n", sep = "")
    print(unique(meta[[col]]))
  }
}
cat("\n")

## ---- 4. EDIT THESE to match your actual metadata columns --------------------
## Set these after reading the diagnostic output above.

SAMPLE_COL     <- "sample"         # after colnames(meta)[1:2] <- c("sample","condition") below
TREATMENT_COL  <- "condition"
ROSI_LEVEL     <- "Rosi"
VEHICLE_LEVEL  <- "DMSO"
PATIENT_COL    <- "patient"        # extracted from sample name below (e.g. "P1DMSO1" -> "P1");
                                    # metadata.csv has no separate patient column, so it must be derived

colnames(meta)[1:2] <- c("sample", "condition")
meta$patient <- sub("^(P[0-9]+).*", "\\1", meta$sample)  # verified pattern against real sample names

if (!SAMPLE_COL %in% colnames(meta) || !TREATMENT_COL %in% colnames(meta)) {
  stop("SAMPLE_COL or TREATMENT_COL not found in metadata. Re-check the printed columns above ",
       "and edit these variables before proceeding.")
}

is_paired <- !is.na(PATIENT_COL) && PATIENT_COL %in% colnames(meta)
cat("Design detected as:", ifelse(is_paired, "PAIRED (patient-matched)", "UNPAIRED/independent"), "\n")
if (!is_paired) {
  cat("If your samples ARE patient-matched, set PATIENT_COL above to the column identifying\n")
  cat("each donor -- a paired test will substantially increase statistical power here.\n\n")
}

rownames(meta) <- meta[[SAMPLE_COL]]
common_samples <- intersect(sample_cols, rownames(meta))
cat(length(common_samples), "of", length(sample_cols), "expression columns matched to metadata.\n")
if (length(common_samples) < length(sample_cols)) {
  cat("WARNING: some samples were not matched -- check for naming mismatches before continuing.\n")
}

## 250 gene symbols are duplicated in this file (confirmed via diagnostic) --
## column_to_rownames() cannot handle duplicate names, and blindly summing
## duplicates would repeat the same error as the earlier microarray dataset.
## Fix: keep every row distinct via make.unique() for the technical row name,
## while preserving the ORIGINAL symbol separately for gene-level matching
## (panel genes were already confirmed not to be among the duplicates).
n_dup <- sum(duplicated(expr_raw[[gene_col]]))
if (n_dup > 0) {
  cat("NOTE:", n_dup, "duplicate gene symbols found in the expression file.\n")
  cat("These rows are kept distinct (not summed) via make.unique(); confirm none of your\n")
  cat("panel genes are affected before trusting downstream panel-level results:\n")
  print(intersect(toupper(expr_raw[[gene_col]][duplicated(expr_raw[[gene_col]]) |
                                                 duplicated(expr_raw[[gene_col]], fromLast = TRUE)]),
                   panel_genes_check <- toupper(c("Adipoq","Cebpa","Fabp4","Plin1","Pnpla3","Scd1",
                                                   "Lpl","Cd36","Slc2a4","Pck1","Ppargc1a","Lep","Retn","Pparg"))))
}
original_symbol <- as.character(expr_raw[[gene_col]])
expr_raw[[gene_col]] <- make.unique(original_symbol)
expr <- expr_raw %>% column_to_rownames(gene_col)
expr <- as.matrix(expr[, common_samples])
mode(expr) <- "numeric"
meta <- meta[common_samples, , drop = FALSE]
meta[[TREATMENT_COL]] <- factor(meta[[TREATMENT_COL]], levels = c(VEHICLE_LEVEL, ROSI_LEVEL))
stopifnot(all(colnames(expr) == rownames(meta)))

group <- meta[[TREATMENT_COL]]
cat("\nGroup assignment:\n"); print(table(group))

## ---- 5. Genome-wide differential expression ----------------------------------
## Branches automatically based on the diagnostic in Section 2. Read the
## printed message to see which branch actually ran.

RUN_LIMMA <- requireNamespace("limma", quietly = TRUE)
DESEQ2_AVAILABLE <- requireNamespace("DESeq2", quietly = TRUE)
looks_like_raw_counts <- all(expr == round(expr), na.rm = TRUE) &&
                          (sd(colSums(expr))/mean(colSums(expr)) > 0.05)

cat("\nlooks_like_raw_counts =", looks_like_raw_counts, "\n")
cat("This determines which differential expression approach is used below.\n")
cat("If this looks wrong given the Section 2 diagnostics, STOP and investigate before trusting results.\n\n")

if (looks_like_raw_counts && DESEQ2_AVAILABLE) {
  # ---- 5A. Raw-count branch: DESeq2 ----
  library(DESeq2)
  design_formula <- if (is_paired) as.formula(paste0("~", PATIENT_COL, " + ", TREATMENT_COL)) else
                                     as.formula(paste0("~", TREATMENT_COL))
  dds <- DESeqDataSetFromMatrix(countData = round(expr), colData = meta, design = design_formula)
  keep <- rowSums(counts(dds) >= 10) >= min(table(group))
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  de_results <- results(dds, contrast = c(TREATMENT_COL, ROSI_LEVEL, VEHICLE_LEVEL)) %>%
    as.data.frame() %>% rownames_to_column("Gene") %>% arrange(padj) %>%
    rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj)
  log_expr <- log2(counts(dds, normalized = TRUE) + 1)
  cat("Used DESeq2 (raw counts detected).\n")

} else if (RUN_LIMMA) {
  # ---- 5B. Pre-normalised branch: limma (as in the Wang et al. pipeline) ----
  library(limma)
  log_expr <- log2(expr + 1)
  if (is_paired) {
    design <- model.matrix(as.formula(paste0("~", PATIENT_COL, " + ", TREATMENT_COL)), data = meta)
  } else {
    design <- model.matrix(~ group)
  }
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit)
  coef_name <- tail(colnames(design), 1)  # the treatment coefficient is last regardless of design
  de_results <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene") %>% rename(t_stat = t)
  cat("Used limma on log2 values (pre-normalised data assumed; DESeq2 not used because",
      "looks_like_raw_counts was FALSE or DESeq2 unavailable).\n")

} else {
  cat("Neither DESeq2 (for raw counts) nor limma is available in this session.\n")
  cat("Falling through to the paired t-test fallback below (Section 5C).\n\n")
}

## ---- 5C. FALLBACK: paired test on patient-level means (only used if neither
##          DESeq2 nor limma is installed; for demonstration/verification only)
## IMPORTANT: with 5 patients x 3 technical replicates x 2 conditions, treating
## all 30 samples as independent would pseudo-replicate within patient. This
## fallback averages the 3 replicates per patient x condition FIRST (giving 5
## paired DMSO/Rosi values per gene), then runs a paired t-test across
## patients -- respecting the true unit of replication. limma on your own
## machine handles this correctly via the patient-blocked design matrix
## without needing this manual averaging step; use these numbers for
## demonstration only, not the manuscript.

if (!exists("de_results")) {
  cat("*** Neither DESeq2 nor limma available in this session. ***\n")
  cat("*** Using patient-averaged PAIRED t-test fallback (demonstration only). ***\n\n")
  log_expr <- log2(expr + 1)

  patient_ids <- unique(meta[[PATIENT_COL]])
  avg_dmso <- sapply(patient_ids, function(p) {
    cols <- rownames(meta)[meta[[PATIENT_COL]] == p & meta[[TREATMENT_COL]] == VEHICLE_LEVEL]
    rowMeans(log_expr[, cols, drop = FALSE])
  })
  avg_rosi <- sapply(patient_ids, function(p) {
    cols <- rownames(meta)[meta[[PATIENT_COL]] == p & meta[[TREATMENT_COL]] == ROSI_LEVEL]
    rowMeans(log_expr[, cols, drop = FALSE])
  })
  colnames(avg_dmso) <- patient_ids; colnames(avg_rosi) <- patient_ids

  de_results <- data.frame(
    Gene = rownames(log_expr),
    logFC = rowMeans(avg_rosi) - rowMeans(avg_dmso),
    AveExpr = rowMeans(cbind(avg_dmso, avg_rosi)),
    stringsAsFactors = FALSE
  )
  pvals <- sapply(seq_len(nrow(log_expr)), function(i) {
    tryCatch(t.test(avg_rosi[i, ], avg_dmso[i, ], paired = TRUE)$p.value, error = function(e) NA)
  })
  de_results$P.Value <- pvals
  de_results$adj.P.Val <- p.adjust(pvals, method = "BH")
  de_results <- de_results %>% arrange(P.Value)
  cat("Fallback used: paired t-test on", length(patient_ids), "patient-level means (n=5 pairs).\n\n")
}

write.csv(de_results, file.path(out_dir, "genomewide_DE_results.csv"), row.names = FALSE)

pnpla3_res <- de_results %>% filter(toupper(Gene) == "PNPLA3")
if (nrow(pnpla3_res) > 0) {
  cat("\n=== PRIMARY RESULT ===\n")
  cat(sprintf("PNPLA3: log2FC = %.3f, adj.P = %.4g (expect POSITIVE logFC: Rosi should INDUCE PNPLA3)\n",
              pnpla3_res$logFC[1], pnpla3_res$adj.P.Val[1]))
} else {
  cat("\n*** PNPLA3 not found by exact symbol match -- check gene ID format (symbol vs Ensembl) ***\n")
}

## ---- 6. Volcano plot (distinct palette) --------------------------------------

volc_df <- de_results %>%
  filter(!is.na(P.Value)) %>%
  mutate(neglog10p = -log10(P.Value),
         sig = case_when(
           adj.P.Val < 0.05 & logFC >  0.5 ~ "Up with Rosi",
           adj.P.Val < 0.05 & logFC < -0.5 ~ "Down with Rosi",
           TRUE ~ "NS"),
         is_pnpla3 = toupper(Gene) == "PNPLA3")

p_volcano <- ggplot(volc_df, aes(x = logFC, y = neglog10p)) +
  geom_point(aes(color = sig), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up with Rosi" = COL_ROSI, "Down with Rosi" = COL_VEHICLE, "NS" = "grey85")) +
  geom_point(data = filter(volc_df, is_pnpla3), color = "black", fill = "#FFD700", shape = 23, size = 4) +
  geom_text_repel(data = filter(volc_df, is_pnpla3), aes(label = "PNPLA3"), fontface = "bold", size = 5) +
  theme_minimal(base_size = 12) +
  labs(title = "Transcriptional response to rosiglitazone in patient ASC-derived adipocytes",
       subtitle = "Hu et al. 2019 (Cell Stem Cell) reanalysis",
       x = "log2 fold change (Rosi / Vehicle)", y = "-log10(P value)", color = NULL)

ggsave(file.path(out_dir, "volcano_ASCRosi.pdf"), p_volcano, width = 7.5, height = 6.5)

## ---- 6b. MA plot, genome-wide, PNPLA3 highlighted (distinct palette) --------
## Confirms PNPLA3's change is not an artefact of low average expression.

p_ma <- ggplot(volc_df, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = sig), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up with Rosi" = COL_ROSI, "Down with Rosi" = COL_VEHICLE, "NS" = "grey85")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(data = filter(volc_df, is_pnpla3), color = "black", fill = "#FFD700", shape = 23, size = 4) +
  geom_text_repel(data = filter(volc_df, is_pnpla3), aes(label = "PNPLA3"), fontface = "bold", size = 5) +
  theme_minimal(base_size = 12) +
  labs(title = "MA plot: PNPLA3 induction is not explained by low average expression",
       x = "Average log2 expression", y = "log2 fold change (Rosi / Vehicle)", color = NULL)

ggsave(file.path(out_dir, "MAplot_ASCRosi.pdf"), p_ma, width = 7.5, height = 6.5)
cat("MA plot written: MAplot_ASCRosi.pdf\n\n")

## ---- 7. PCA (distinct palette) ------------------------------------------------

var_genes <- order(apply(log_expr, 1, var), decreasing = TRUE)[1:min(2000, nrow(log_expr))]
pca <- prcomp(t(log_expr[var_genes, ]), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>% rownames_to_column("sample") %>% left_join(
  meta, by = "sample")
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[TREATMENT_COL]])) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL))) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA: rosiglitazone effect across patient-derived adipocyte lines",
       x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"), color = NULL)
if (is_paired) {
  p_pca <- p_pca + geom_line(aes(group = .data[[PATIENT_COL]]), color = "grey70", alpha = 0.5)
}
ggsave(file.path(out_dir, "PCA_ASCRosi.pdf"), p_pca, width = 7, height = 6)

## ---- 8. Curated adipocyte-program panel (human symbols; PPARG included) -----
## Unlike the Wang et al. mouse microarray (where 3 conflicting Pparg probes
## required exclusion), this RNA-seq matrix has exactly one row for PPARG and
## one for PNPLA3 (confirmed: neither appears among this file's ~250 duplicate
## gene symbols), so PPARG is included here with no ambiguity to resolve.

panel_genes <- toupper(c("Adipoq","Cebpa","Fabp4","Plin1","Pnpla3","Scd1","Lpl","Cd36",
                          "Slc2a4","Pck1","Ppargc1a","Lep","Retn","Pparg"))

panel_de <- de_results %>% filter(toupper(Gene) %in% panel_genes)
write.csv(panel_de, file.path(out_dir, "adipocyte_panel_DE_results.csv"), row.names = FALSE)
cat("\n=== Adipocyte panel DE results ===\n")
print(panel_de %>% select(Gene, logFC, P.Value, adj.P.Val))

## ---- 9. Forest plot (distinct style: diamonds, teal/coral) -------------------

df_forest <- panel_de %>% arrange(logFC) %>%
  mutate(is_pnpla3 = toupper(Gene) == "PNPLA3")

p_forest <- ggplot(df_forest, aes(x = logFC, y = reorder(Gene, logFC))) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(color = is_pnpla3, shape = is_pnpla3), size = 4) +
  scale_color_manual(values = c("TRUE" = "#FFB300", "FALSE" = "#2C7FB8"), guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  theme_minimal(base_size = 11) +
  labs(title = "Effect of rosiglitazone on the adipocyte gene panel",
       x = "log2 fold change (Rosi vs Vehicle)", y = NULL)

ggsave(file.path(out_dir, "forest_plot_ASCRosi.pdf"), p_forest, width = 7, height = 6)

## ---- 10. PNPLA3 dot plot (violin + jitter, distinct from mouse boxplots) -----

pnpla3_gene <- de_results %>% filter(toupper(Gene) == "PNPLA3") %>% pull(Gene) %>% first()
pnpla3_vals <- data.frame(sample = colnames(log_expr), expr = log_expr[pnpla3_gene, ]) %>%
  left_join(meta, by = "sample")

pnpla3_stats <- pnpla3_vals %>% group_by(.data[[TREATMENT_COL]]) %>%
  summarise(mean = mean(expr), sem = sd(expr)/sqrt(n()), n = n())

pnpla3_pval <- de_results %>% filter(toupper(Gene) == "PNPLA3") %>% pull(P.Value) %>% first()

p_pnpla3 <- ggplot(pnpla3_vals, aes(x = .data[[TREATMENT_COL]], y = expr, fill = .data[[TREATMENT_COL]])) +
  geom_violin(alpha = 0.5, color = NA) +
  geom_jitter(width = 0.08, size = 2, color = "grey20", alpha = 0.7) +
  geom_pointrange(data = pnpla3_stats, inherit.aes = FALSE,
                   aes(x = .data[[TREATMENT_COL]], y = mean, ymin = mean - sem, ymax = mean + sem),
                   color = "black", size = 0.9, fatten = 2) +
  scale_fill_manual(values = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL))) +
  theme_minimal(base_size = 13) +
  labs(title = "PNPLA3 expression: Vehicle vs Rosiglitazone",
       subtitle = paste0("n = ", nrow(pnpla3_vals), " samples (",
                          ifelse(is_paired, paste0(length(unique(meta[[PATIENT_COL]])), " patients, paired design"),
                                 "independent samples"),
                          "); P = ", signif(pnpla3_pval, 3),
                          " (see Section 5 for which test/branch was used); black points = mean \u00B1 SEM"),
       x = NULL, y = "log2 expression (PNPLA3)")
if (is_paired) {
  p_pnpla3 <- p_pnpla3 + geom_line(aes(group = .data[[PATIENT_COL]]), color = "grey60", alpha = 0.4)
}
ggsave(file.path(out_dir, "PNPLA3_violinplot.pdf"), p_pnpla3, width = 5.5, height = 6)

## ---- 10b. Individual violin plots for EVERY panel gene (one file each) ------
## Same visual style as the PNPLA3 plot above (violin + jitter + mean/-SEM +
## paired lines), generalised into a function and applied to the whole
## adipocyte panel. Each gene gets its own standalone PDF, not a facet grid.

make_gene_violin <- function(gene_upper) {
  row_name <- rownames(expr)[which(toupper(original_symbol) == gene_upper)[1]]
  if (is.na(row_name)) {
    cat("Skipping", gene_upper, "-- not found by exact symbol match.\n")
    return(invisible(NULL))
  }

  gene_vals <- data.frame(sample = colnames(log_expr), expr = as.numeric(log_expr[row_name, ])) %>%
    left_join(meta, by = "sample")

  gene_stats <- gene_vals %>% group_by(.data[[TREATMENT_COL]]) %>%
    summarise(mean = mean(expr), sem = sd(expr)/sqrt(n()), n = n())

  gene_pval <- de_results %>% filter(toupper(Gene) == gene_upper) %>% pull(P.Value) %>% first()
  pval_label <- if (length(gene_pval) == 0 || is.na(gene_pval)) "NA" else signif(gene_pval, 3)

  p <- ggplot(gene_vals, aes(x = .data[[TREATMENT_COL]], y = expr, fill = .data[[TREATMENT_COL]])) +
    geom_violin(alpha = 0.5, color = NA) +
    geom_jitter(width = 0.08, size = 2, color = "grey20", alpha = 0.7) +
    geom_pointrange(data = gene_stats, inherit.aes = FALSE,
                     aes(x = .data[[TREATMENT_COL]], y = mean, ymin = mean - sem, ymax = mean + sem),
                     color = "black", size = 0.9, fatten = 2) +
    scale_fill_manual(values = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL))) +
    theme_minimal(base_size = 13) +
    labs(title = paste0(gene_upper, " expression: Vehicle vs Rosiglitazone"),
         subtitle = paste0("n = ", nrow(gene_vals), " samples (",
                            ifelse(is_paired, paste0(length(unique(meta[[PATIENT_COL]])), " patients, paired design"),
                                   "independent samples"),
                            "); P = ", pval_label, "; black points = mean \u00B1 SEM"),
         x = NULL, y = paste0("log2 expression (", gene_upper, ")"))
  if (is_paired) {
    p <- p + geom_line(aes(group = .data[[PATIENT_COL]]), color = "grey60", alpha = 0.4)
  }
  ggsave(file.path(out_dir, paste0(gene_upper, "_violinplot.pdf")), p, width = 6, height = 6.5)
  cat("Violin plot written:", paste0(gene_upper, "_violinplot.pdf"), " (P =", pval_label, ")\n")
}

cat("\n=== Generating individual violin plots for every panel gene ===\n")
invisible(lapply(panel_genes, make_gene_violin))
cat("\n")

## ---- 11. Heatmap (distinct color scale: blue-white-coral) --------------------

panel_present <- intersect(panel_genes, toupper(rownames(log_expr)))
heat_rows <- rownames(log_expr)[toupper(rownames(log_expr)) %in% panel_present]
heat_mat <- log_expr[heat_rows, ]
heat_z <- t(scale(t(heat_mat)))

col_anno <- data.frame(Treatment = meta[[TREATMENT_COL]]); rownames(col_anno) <- rownames(meta)
anno_colors <- list(Treatment = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL)))

pdf(file.path(out_dir, "adipocyte_program_heatmap_ASCRosi.pdf"), width = 7, height = 8)
pheatmap(heat_z, cluster_rows = TRUE, cluster_cols = is_paired == FALSE,
         annotation_col = col_anno, annotation_colors = anno_colors,
         color = colorRampPalette(PALETTE_DIV)(100),
         main = "Adipocyte transcriptional program (Z-score)\nRosiglitazone vs Vehicle, patient ASC-adipocytes",
         fontsize_row = 9, border_color = NA)
dev.off()

## ---- 12. Composite adipocyte-program score (PNPLA3 excluded from score) -----

score_genes <- setdiff(panel_present, "PNPLA3")
score_rows <- rownames(log_expr)[toupper(rownames(log_expr)) %in% score_genes]
score_z <- t(scale(t(log_expr[score_rows, ])))
program_score <- colMeans(score_z, na.rm = TRUE)

score_df <- data.frame(sample = names(program_score), program_score = program_score) %>%
  left_join(meta, by = "sample") %>%
  left_join(pnpla3_vals %>% select(sample, pnpla3_expr = expr), by = "sample")
write.csv(score_df, file.path(out_dir, "program_score_ASCRosi.csv"), row.names = FALSE)

p_score <- ggplot(score_df, aes(x = program_score, y = pnpla3_expr, color = .data[[TREATMENT_COL]])) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = program_score, y = pnpla3_expr),
              se = TRUE, color = "grey30", linetype = "dashed") +
  scale_color_manual(values = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL))) +
  theme_minimal(base_size = 12) +
  labs(title = "PNPLA3 expression vs adipocyte-program score",
       subtitle = "Covariation with the Rosi-driven adipocyte program, not evidence of direct regulation",
       x = "Adipocyte program score (Z, panel excl. PNPLA3)", y = "PNPLA3 (log2 expression)", color = NULL)
ggsave(file.path(out_dir, "PNPLA3_vs_program_score_ASCRosi.pdf"), p_score, width = 7, height = 5.5)

## ---- 12b. Ranked log2FC bar plot: PNPLA3 relative to the panel --------------
## Distinct from the forest plot (which shows CIs) -- this is the simple
## magnitude comparison recommended, PNPLA3 highlighted.

p_ranked <- ggplot(df_forest, aes(x = reorder(Gene, logFC), y = logFC)) +
  geom_col(aes(fill = is_pnpla3)) +
  scale_fill_manual(values = c("TRUE" = "#FFB300", "FALSE" = "grey65"), guide = "none") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "Magnitude of PNPLA3 induction by rosiglitazone,\nrelative to established adipocyte genes",
       x = NULL, y = "log2 fold change (Rosi vs Vehicle)")
ggsave(file.path(out_dir, "ranked_log2FC_ASCRosi.pdf"), p_ranked, width = 7, height = 6)
cat("Ranked log2FC plot written: ranked_log2FC_ASCRosi.pdf\n\n")

## ---- 12c. Pairwise PNPLA3 correlation scatterplots (NOT a network) ----------
## All 30 individual samples shown, coloured by treatment. NOTE: because
## samples are patient-matched (3 technical replicates x 5 patients x 2
## conditions), these 30 points are not 30 independent biological units --
## much of any correlation reflects the shared Rosi/DMSO split across
## patients. This is shown as a visual of covariation with the coordinated
## adipocyte response, not as evidence that any one gene regulates another.

partner_genes <- intersect(toupper(c("Adipoq","Fabp4","Cd36","Slc2a4","Lpl")),
                            toupper(original_symbol))
pnpla3_row_name <- rownames(expr)[which(toupper(original_symbol) == "PNPLA3")[1]]

scatter_df <- lapply(partner_genes, function(g) {
  row_name <- rownames(expr)[which(toupper(original_symbol) == g)[1]]
  data.frame(sample = colnames(log_expr),
             partner_gene = g,
             partner_expr = as.numeric(log_expr[row_name, ]),
             pnpla3_expr  = as.numeric(log_expr[pnpla3_row_name, ]))
}) %>% bind_rows() %>% left_join(meta, by = "sample")

cor_stats <- scatter_df %>% group_by(partner_gene) %>%
  summarise(r = cor(pnpla3_expr, partner_expr, method = "pearson"),
            p = cor.test(pnpla3_expr, partner_expr)$p.value, n = n())
cat("=== PNPLA3 vs individual adipocyte genes (n=30 samples, patient-matched -- interpret as covariation) ===\n")
print(cor_stats)
write.csv(cor_stats, file.path(out_dir, "pnpla3_pairwise_correlations_ASCRosi.csv"), row.names = FALSE)

p_scatter <- ggplot(scatter_df, aes(x = partner_expr, y = pnpla3_expr, color = .data[[TREATMENT_COL]])) +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_color_manual(values = setNames(c(COL_VEHICLE, COL_ROSI), c(VEHICLE_LEVEL, ROSI_LEVEL))) +
  facet_wrap(~partner_gene, scales = "free_x") +
  theme_minimal(base_size = 10) +
  labs(title = "PNPLA3 vs individual adipocyte genes across all 30 samples",
       subtitle = "Patient-matched design -- shown as covariation with the coordinated Rosi response, not direct regulation",
       x = "Partner gene (log2 expression)", y = "PNPLA3 (log2 expression)", color = NULL)
ggsave(file.path(out_dir, "PNPLA3_pairwise_scatterplots_ASCRosi.pdf"), p_scatter, width = 9, height = 6)
cat("Pairwise scatterplots written: PNPLA3_pairwise_scatterplots_ASCRosi.pdf\n\n")

## ---- 12d. Descriptive tiers (narrative aid only, not a separate test) -------

tiers <- panel_de %>%
  mutate(tier = case_when(
    adj.P.Val < 0.05 & abs(logFC) > 1   ~ "Strongly altered",
    adj.P.Val < 0.05 & abs(logFC) <= 1  ~ "Moderately altered",
    TRUE                                 ~ "Relatively preserved"
  )) %>%
  select(Gene, logFC, adj.P.Val, tier) %>% arrange(tier, desc(logFC))

write.csv(tiers, file.path(out_dir, "descriptive_tiers_ASCRosi.csv"), row.names = FALSE)
cat("=== Descriptive tiers (narrative aid only) ===\n")
print(tiers)
cat("\nCaution: genes in 'Relatively preserved' were not strongly altered under THIS specific\n")
cat("treatment/timepoint/context -- this does not mean rosiglitazone/PPARG never regulates them.\n\n")

## ---- 13. GSEA -- Hallmark, HUMAN native gene sets -----------------------------
## Uses lessons learned from the mouse pipeline: explicit collection code and
## db_species, checked against msigdbr_collections() rather than assumed.

if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE) &&
    "t_stat" %in% colnames(de_results)) {
  library(fgsea); library(msigdbr)

  cat("\nChecking available human collections before querying (see prior mouse-pipeline lesson:\n")
  cat("mouse required collection='MH' for db_species='MM', NOT the naive 'H') ...\n")
  avail_collections <- msigdbr_collections(db_species = "HS")
  print(avail_collections)
  if (!"H" %in% avail_collections$gs_collection) {
    stop("Collection 'H' not found for db_species='HS' in this msigdbr version. ",
         "Check the printed table above for the correct Hallmark identifier and update the ",
         "collection= argument below accordingly before proceeding.")
  }
  cat("Confirmed: 'H' is a valid collection for db_species='HS'.\n\n")

  ranked <- de_results %>%
    filter(!is.na(Gene), Gene != "") %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(t_stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>% arrange(desc(t_stat))
  rank_vec <- setNames(ranked$t_stat, toupper(ranked$Gene))

  hallmark <- msigdbr(species = "Homo sapiens", collection = "H", db_species = "HS")
  pathways_hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)

  gsea_res <- fgsea(pathways = pathways_hallmark, stats = rank_vec, minSize = 10, maxSize = 500) %>%
    arrange(padj)
  write.csv(as.data.frame(gsea_res %>% select(-leadingEdge)),
            file.path(out_dir, "GSEA_hallmark_ASCRosi.csv"), row.names = FALSE)

  relevant <- gsea_res %>%
    filter(grepl("ADIPOGENESIS|FATTY_ACID|PPAR|OXIDATIVE_PHOSPHORYLATION|GLYCOLYSIS", pathway, ignore.case = TRUE))
  p_gsea <- ggplot(relevant, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE" = COL_ROSI, "FALSE" = COL_VEHICLE), guide = "none") +
    theme_minimal(base_size = 11) +
    labs(title = "GSEA: adipocyte/metabolic pathways induced by rosiglitazone",
         x = NULL, y = "Normalized Enrichment Score (NES)")
  ggsave(file.path(out_dir, "GSEA_barplot_ASCRosi.pdf"), p_gsea, width = 9, height = 5)

} else {
  cat("\n*** GSEA skipped: requires fgsea, msigdbr, AND a limma-derived 't_stat' column ***\n")
  cat("*** (i.e. the DESeq2 branch was used, or these packages are missing).        ***\n")
  cat("*** If DESeq2 was used, rank by the 'stat' column instead and adapt this     ***\n")
  cat("*** section's group_by/slice_max call accordingly before running GSEA.       ***\n\n")
}

## ---- 14. Summary ---------------------------------------------------------------

cat("\n\n===== DONE (UNTESTED SCRIPT -- verify every diagnostic message above) =====\n")
cat("Outputs written to:", out_dir, "\n")
cat("Before trusting any result: re-check Sections 2-4 diagnostics matched your expectations,\n")
cat("confirm PATIENT_COL/paired-design detection is correct, and confirm the looks_like_raw_counts\n")
cat("branch decision in Section 5 was appropriate for your actual data.\n")
