# --- Portable paths: set these two lines (or run from the project root with data/ and results/ folders) ---
DATA_DIR <- Sys.getenv("PPARG_DATA", "data"); RESULTS_DIR <- Sys.getenv("PPARG_RESULTS", "results")

#!/usr/bin/env Rscript
# ==============================================================================
# PPARg deletion -> Pnpla3 and the adipocyte transcriptional program
# Dataset: GSE62937 (Wang et al. 2015, Nat Cell Biol), GPL6887 (Illumina
#          MouseWG-6 v2.0 BeadChip). Adn-Pparg-KO vs Control, n=3 vs n=3
#          (each sample pooled from 9-14 mice per the paper's Methods).
#
# WHAT THIS SCRIPT DOES AND WHY, BEFORE YOU RUN ANYTHING:
#
#  - This is genome-wide, probe-level data extracted directly from the GEO
#    series matrix (45,281 probes), annotated by gene symbol using the
#    official GPL6887 manifest. This is NOT the same as the earlier
#    "countmatrix.csv" gene-symbol-collapsed file -- this script starts from
#    probe IDs and never collapses/sums multi-probe genes.
#
#  - No raw CEL files are available, so RMA normalisation from scratch is not
#    possible. The values here are the already-normalised intensities
#    deposited in GEO (per the paper's Methods: GeneSpring GX normalisation).
#    limma is used directly on these normalised, log2-transformed values,
#    which is standard practice when only normalised intensities are
#    available. This does NOT re-derive normalisation; it borrows limma's
#    moderated statistics (empirical Bayes variance shrinkage), which is
#    exactly what benefits a small n=3 vs n=3 design.
#
#  - limma, fgsea and msigdbr are Bioconductor/CRAN packages that could not be
#    installed in the sandbox used to build/verify this script (no internet
#    access to bioconductor.org or cran.r-project.org in that environment).
#    Sections 2 (genome-wide DE) and 9 (GSEA) THEREFORE MUST BE RUN ON YOUR
#    OWN MACHINE, where install.packages()/BiocManager::install() will work
#    normally. All other sections (PCA, volcano, MA plot, forest plot, Pnpla3
#    dot plot, heatmap, composite score, correlation scatterplots, ranked
#    plot) use only base R + CRAN packages available via apt, and were run
#    and verified against your real data before this script was delivered.
#
#  - PPARG ITSELF IS DELIBERATELY EXCLUDED from the curated adipocyte panel
#    throughout (heatmap, forest plot, composite score, etc.), consistent
#    with prior analysis: 3 Illumina probes map to Pparg and disagree in
#    direction (2 up, 1 down, in the labelled KO), for reasons traced to the
#    probe level but not fully resolved (see project notes). Do not add Pparg
#    back into these panels without addressing that first.
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

## ---- 1. Paths (EDIT IF NEEDED) ----------------------------------------------

series_matrix_path <- file.path(DATA_DIR, "GSE62937_series_matrix.txt")
manifest_path       <- file.path(DATA_DIR, "GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt")
metadata_path       <- file.path(DATA_DIR, "metadata.csv")
out_dir             <- file.path(RESULTS_DIR, "mouse_genomewide")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_SAMPLES <- c("GSM1536683","GSM1536684","GSM1536685",
                    "GSM1536686","GSM1536687","GSM1536688")

## ---- 2. Load genome-wide series matrix, extract our 6 samples ---------------

sm_lines <- readLines(series_matrix_path)
header_line <- grep("^!Sample_geo_accession", sm_lines, value = TRUE)
header <- strsplit(gsub('"', '', header_line), "\t")[[1]]

data_lines <- grep("^ILMN_", sm_lines, value = TRUE)
mat_raw <- read.delim(text = paste(data_lines, collapse = "\n"), header = FALSE,
                       stringsAsFactors = FALSE)
colnames(mat_raw) <- header
rownames(mat_raw) <- mat_raw[[1]]

expr <- as.matrix(mat_raw[, TARGET_SAMPLES])
mode(expr) <- "numeric"
cat("Genome-wide matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

## ---- 3. Load manifest, annotate probes with gene symbols --------------------

manifest <- read.delim(manifest_path, skip = 8, stringsAsFactors = FALSE)
manifest <- manifest %>% select(Probe_Id, Symbol, RefSeq_ID) %>% filter(!is.na(Probe_Id), Probe_Id != "")
rownames(manifest) <- manifest$Probe_Id

gene_symbol <- manifest[rownames(expr), "Symbol"]
cat("Probes annotated with a gene symbol:", sum(!is.na(gene_symbol)), "of", length(gene_symbol), "\n\n")

## ---- 4. Metadata / design ----------------------------------------------------

meta <- read.csv(metadata_path, check.names = FALSE)
colnames(meta)[1:2] <- c("sample", "condition")
meta$condition <- factor(meta$condition, levels = c("Con Adn-PPARG-/-", "Adn-PPARG-/-"))
meta <- meta[match(TARGET_SAMPLES, meta$sample), ]
stopifnot(all(colnames(expr) == meta$sample))

group <- meta$condition
cat("Group assignment:\n"); print(table(group))

log_expr <- log2(expr + 1)  # values are already normalised intensities; log2 for variance stabilisation

# ==============================================================================
## ---- 5. GENOME-WIDE DIFFERENTIAL EXPRESSION (limma) -- RUN ON YOUR OWN MACHINE
## This section requires Bioconductor's limma and could not be executed in the
## sandbox used to build this script. Install locally with:
##   if (!require("BiocManager")) install.packages("BiocManager")
##   BiocManager::install("limma")
## Then run this block as-is.
# ==============================================================================

RUN_LIMMA <- requireNamespace("limma", quietly = TRUE)

if (RUN_LIMMA) {
  library(limma)
  design <- model.matrix(~ group)
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit)
  limma_results <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
    rownames_to_column("Probe_Id") %>%
    rename(t_stat = t) %>%
    mutate(Symbol = manifest[Probe_Id, "Symbol"]) %>%
    relocate(Symbol, .after = Probe_Id)
  write.csv(limma_results, file.path(out_dir, "genomewide_limma_results.csv"), row.names = FALSE)
  cat("limma genome-wide results written: genomewide_limma_results.csv\n")

  pnpla3_limma <- limma_results %>% filter(Symbol == "Pnpla3")
  cat("\n=== PRIMARY RESULT (limma, moderated statistics) ===\n")
  cat(sprintf("Pnpla3: log2FC = %.3f (fold change = %.2fx), adj.P = %.4g\n",
              pnpla3_limma$logFC, 2^pnpla3_limma$logFC, pnpla3_limma$adj.P.Val))

} else {
  cat("\n*** limma is not installed in this R session. ***\n")
  cat("*** Install it (BiocManager::install('limma')) and re-run this script ***\n")
  cat("*** to get the genome-wide differential expression table, volcano/MA  ***\n")
  cat("*** plots based on moderated statistics, and the GSEA ranking below.  ***\n\n")
  cat("Falling back to per-probe Welch's t-test genome-wide for demonstration\n")
  cat("purposes ONLY -- these are NOT moderated statistics and should NOT be\n")
  cat("used as the manuscript's primary result. Re-run with limma installed.\n\n")

  # Fallback: plain per-probe Welch t-test (for demonstration/plotting only)
  ctrl_idx <- which(group == "Con Adn-PPARG-/-")
  ko_idx   <- which(group == "Adn-PPARG-/-")

  fallback_results <- data.frame(
    Probe_Id = rownames(log_expr),
    Symbol   = gene_symbol,
    logFC    = rowMeans(log_expr[, ko_idx]) - rowMeans(log_expr[, ctrl_idx]),
    AveExpr  = rowMeans(log_expr),
    stringsAsFactors = FALSE
  )
  pvals <- apply(log_expr, 1, function(x) {
    tryCatch(t.test(x[ko_idx], x[ctrl_idx], var.equal = FALSE)$p.value, error = function(e) NA)
  })
  fallback_results$P.Value <- pvals
  fallback_results$adj.P.Val <- p.adjust(pvals, method = "BH")
  fallback_results <- fallback_results %>% arrange(P.Value)
  write.csv(fallback_results, file.path(out_dir, "genomewide_WELCH_FALLBACK_results.csv"), row.names = FALSE)
  limma_results <- fallback_results  # so downstream sections still run for demonstration

  pnpla3_fb <- fallback_results %>% filter(Symbol == "Pnpla3")
  cat("=== FALLBACK RESULT (Welch t-test, NOT for manuscript use) ===\n")
  cat(sprintf("Pnpla3: log2FC = %.3f (fold change = %.2fx), adj.P = %.4g\n",
              pnpla3_fb$logFC, 2^pnpla3_fb$logFC, pnpla3_fb$adj.P.Val))
}

## ---- 6. Volcano plot, genome-wide, Pnpla3 highlighted -----------------------

volc_df <- limma_results %>%
  mutate(neglog10p = -log10(P.Value),
         sig = case_when(
           adj.P.Val < 0.05 & logFC >  0.5 ~ "Up in KO",
           adj.P.Val < 0.05 & logFC < -0.5 ~ "Down in KO",
           TRUE ~ "NS"),
         is_pnpla3 = Symbol == "Pnpla3")

p_volcano <- ggplot(volc_df, aes(x = logFC, y = neglog10p)) +
  geom_point(aes(color = sig), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up in KO" = "firebrick", "Down in KO" = "steelblue", "NS" = "grey80")) +
  geom_point(data = filter(volc_df, is_pnpla3), color = "black", fill = "gold", shape = 21, size = 3.5) +
  geom_text_repel(data = filter(volc_df, is_pnpla3), aes(label = "Pnpla3"), fontface = "bold", size = 5) +
  theme_bw(base_size = 12) +
  labs(title = "Genome-wide transcriptional response to Pparg deletion (GSE62937)",
       x = "log2 fold change (Adn-Pparg-KO / Control)", y = "-log10(P value)", color = NULL)

ggsave(file.path(out_dir, "volcano_genomewide.pdf"), p_volcano, width = 7.5, height = 6.5)
cat("\nVolcano plot written: volcano_genomewide.pdf\n")

## ---- 7. MA plot, genome-wide, Pnpla3 highlighted ----------------------------

p_ma <- ggplot(volc_df, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = sig), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up in KO" = "firebrick", "Down in KO" = "steelblue", "NS" = "grey80")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = filter(volc_df, is_pnpla3), color = "black", fill = "gold", shape = 21, size = 3.5) +
  geom_text_repel(data = filter(volc_df, is_pnpla3), aes(label = "Pnpla3"), fontface = "bold", size = 5) +
  theme_bw(base_size = 12) +
  labs(title = "MA plot: Pnpla3 change is not explained by low average expression",
       x = "Average log2 expression", y = "log2 fold change (KO / Control)", color = NULL)

ggsave(file.path(out_dir, "MAplot_genomewide.pdf"), p_ma, width = 7.5, height = 6.5)
cat("MA plot written: MAplot_genomewide.pdf\n\n")

## ---- 8. PCA, genome-wide -----------------------------------------------------
## Uses the most variable probes (standard practice) rather than all 45k,
## which are dominated by noise for lowly-expressed/invariant probes.

var_probes <- order(apply(log_expr, 1, var), decreasing = TRUE)[1:min(2000, nrow(log_expr))]
pca <- prcomp(t(log_expr[var_probes, ]), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(title = "PCA: Adn-Pparg-KO produces a distinct global transcriptional phenotype",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"), color = NULL)

ggsave(file.path(out_dir, "PCA_genomewide.pdf"), p_pca, width = 7, height = 6)
cat("PCA plot written: PCA_genomewide.pdf\n\n")

# ==============================================================================
## ---- 9. GSEA -- RUN ON YOUR OWN MACHINE
## Requires fgsea (Bioconductor) and msigdbr (CRAN); not installable in the
## sandbox used to build this script. Install locally with:
##   BiocManager::install("fgsea")
##   install.packages("msigdbr")
## Ranks every probe by the limma moderated t-statistic (NOT logFC alone --
## t-statistic accounts for variance and is the standard GSEA ranking metric),
## collapses to gene symbol by taking the probe with the largest |t| per gene
## (a common, defensible approach for multi-probe array-based GSEA), and tests
## MSigDB Hallmark + GO Biological Process gene sets.
# ==============================================================================

if (RUN_LIMMA && requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
  library(fgsea); library(msigdbr)

  ranked <- limma_results %>%
    filter(!is.na(Symbol), Symbol != "") %>%
    group_by(Symbol) %>%
    slice_max(order_by = abs(t_stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(t_stat))

  rank_vec <- setNames(ranked$t_stat, ranked$Symbol)

  hallmark <- msigdbr(species = "Mus musculus", category = "H")
  pathways_hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)

  gsea_res <- fgsea(pathways = pathways_hallmark, stats = rank_vec, minSize = 10, maxSize = 500)
  gsea_res <- gsea_res %>% arrange(padj)
  fwrite <- function(x, f) write.csv(as.data.frame(x %>% select(-leadingEdge)), f, row.names = FALSE)
  fwrite(gsea_res, file.path(out_dir, "GSEA_hallmark_results.csv"))
  cat("GSEA (Hallmark) written: GSEA_hallmark_results.csv\n")

  # Highlight the pathways relevant to the biological question
  relevant <- gsea_res %>%
    filter(grepl("ADIPOGENESIS|FATTY_ACID|PPAR|MTORC|OXIDATIVE_PHOSPHORYLATION|GLYCOLYSIS", pathway, ignore.case = TRUE))
  cat("\nRelevant Hallmark pathways:\n")
  print(relevant %>% select(pathway, NES, pval, padj))

  p_gsea <- ggplot(relevant, aes(x = reorder(pathway, NES), y = NES, fill = NES < 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("firebrick", "steelblue"), guide = "none") +
    theme_bw(base_size = 11) +
    labs(title = "GSEA: adipocyte/metabolic Hallmark pathways after Pparg deletion",
         x = NULL, y = "Normalized Enrichment Score (NES)")
  ggsave(file.path(out_dir, "GSEA_barplot.pdf"), p_gsea, width = 9, height = 5)
  cat("GSEA bar plot written: GSEA_barplot.pdf\n\n")

} else if (!RUN_LIMMA) {
  cat("*** GSEA requires limma's moderated t-statistic for ranking. ***\n")
  cat("*** Install limma (see Section 5) and re-run the whole script. ***\n\n")
} else {
  cat("*** fgsea and/or msigdbr not installed. ***\n")
  cat("*** Install (see comments above) and re-run this script for GSEA. ***\n\n")
}

## ---- 10. Curated adipocyte-program panel (Pparg excluded; see header) -------

panel_genes <- c("Adipoq","Cebpa","Fabp4","Plin1",          # adipocyte identity
                  "Pnpla3","Scd1","Lpl","Cd36",              # lipid metabolism
                  "Slc2a4","Pck1","Ppargc1a",                # glucose/metabolic
                  "Lep","Retn")                              # adipokines

panel_de <- limma_results %>% filter(Symbol %in% panel_genes) %>%
  mutate(n_probes_for_gene = ave(Symbol, Symbol, FUN = length))
write.csv(panel_de, file.path(out_dir, "adipocyte_panel_DE_results.csv"), row.names = FALSE)
cat("=== Adipocyte panel DE results (Pparg excluded) ===\n")
print(panel_de %>% select(Symbol, logFC, P.Value, adj.P.Val, n_probes_for_gene))

## ---- 11. Forest / dot plot of log2FC with 95% CI ----------------------------
## CI derived from limma's t and logFC (SE = logFC / t); if running the Welch
## fallback, CI uses the t-test's own CI directly instead.

if (RUN_LIMMA) {
  df_forest <- fit %>% topTable(coef = 2, number = Inf) %>%
    rownames_to_column("Probe_Id") %>%
    rename(t_stat = t) %>%
    mutate(Symbol = manifest[Probe_Id, "Symbol"]) %>%
    filter(Symbol %in% panel_genes) %>%
    mutate(SE = logFC / t_stat,
           CI_lower = logFC - qt(0.975, df = fit$df.total[1]) * SE,
           CI_upper = logFC + qt(0.975, df = fit$df.total[1]) * SE,
           row_label = paste0(Symbol, ifelse(ave(Symbol, Symbol, FUN=length) > 1,
                                              paste0(" (", Probe_Id, ")"), "")))
} else {
  # Fallback CI via per-probe t.test on log-expression directly
  df_forest <- panel_de
  ctrl_idx <- which(group == "Con Adn-PPARG-/-"); ko_idx <- which(group == "Adn-PPARG-/-")
  ci_list <- lapply(df_forest$Probe_Id, function(pid) {
    tt <- t.test(log_expr[pid, ko_idx], log_expr[pid, ctrl_idx], var.equal = FALSE)
    data.frame(CI_lower = tt$conf.int[1], CI_upper = tt$conf.int[2])
  })
  df_forest <- bind_cols(df_forest, bind_rows(ci_list))
  df_forest$row_label <- paste0(df_forest$Symbol,
                                 ifelse(ave(df_forest$Symbol, df_forest$Symbol, FUN=length) > 1,
                                        paste0(" (", df_forest$Probe_Id, ")"), ""))
}

df_forest <- df_forest %>% arrange(logFC)

p_forest <- ggplot(df_forest, aes(x = logFC, y = reorder(row_label, logFC))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_point(aes(color = Symbol == "Pnpla3"), size = 3) +
  scale_color_manual(values = c("TRUE" = "gold3", "FALSE" = "black"), guide = "none") +
  theme_bw(base_size = 11) +
  labs(title = "Effect of Pparg deletion on adipocyte genes",
       x = "log2 fold change (KO vs Control), 95% CI", y = NULL)

ggsave(file.path(out_dir, "forest_plot_panel.pdf"), p_forest, width = 7, height = 6)
cat("\nForest plot written: forest_plot_panel.pdf\n\n")

## ---- 12. Individual violin plots for EVERY panel gene (one file each) -------
## New, third distinct colour identity for this dataset: deep indigo (Control)
## vs warm gold/mustard (KO) -- visually distinct from both the ASC-Rosi
## (teal/coral) and earlier mouse heatmap (blue/orange) palettes.

COL_CONTROL_MOUSE <- "#5B4B8A"   # deep indigo/plum
COL_KO_MOUSE      <- "#E4A62B"   # warm gold/mustard

# Kept in scope for downstream sections (composite score, ranked plot) that
# reference Pnpla3 specifically -- the plot itself is now produced by the
# generalised loop below, which covers Pnpla3 along with every other gene.
pnpla3_probe <- manifest %>% filter(Symbol == "Pnpla3") %>% pull(Probe_Id)
pnpla3_vals <- data.frame(sample = colnames(expr), expr = as.numeric(expr[pnpla3_probe, ])) %>%
  left_join(meta, by = "sample")

make_mouse_gene_violin <- function(gene_symbol_upper_safe) {
  probe_id <- manifest %>% filter(Symbol == gene_symbol_upper_safe) %>% pull(Probe_Id)
  if (length(probe_id) == 0) {
    cat("Skipping", gene_symbol_upper_safe, "-- not found by exact symbol match.\n")
    return(invisible(NULL))
  }
  if (length(probe_id) > 1) {
    cat("Skipping", gene_symbol_upper_safe, "in this loop -- multiple probes, handled separately (see network section).\n")
    return(invisible(NULL))
  }

  gene_vals <- data.frame(sample = colnames(expr), expr = as.numeric(expr[probe_id, ])) %>%
    left_join(meta, by = "sample")

  gene_stats <- gene_vals %>% group_by(condition) %>%
    summarise(mean = mean(expr), sem = sd(expr)/sqrt(n()), n = n())

  gene_pval <- limma_results %>% filter(Symbol == gene_symbol_upper_safe) %>% pull(adj.P.Val) %>% first()
  pval_label <- if (length(gene_pval) == 0 || is.na(gene_pval)) "NA" else paste0("adj.P = ", signif(gene_pval, 3))

  p <- ggplot(gene_vals, aes(x = condition, y = expr, fill = condition)) +
    geom_violin(alpha = 0.55, color = NA, trim = FALSE) +
    geom_jitter(width = 0.06, size = 2.2, color = "grey20", alpha = 0.8) +
    geom_pointrange(data = gene_stats, inherit.aes = FALSE,
                     aes(x = condition, y = mean, ymin = mean - sem, ymax = mean + sem),
                     color = "black", size = 0.9, fatten = 2) +
    scale_fill_manual(values = c("Con Adn-PPARG-/-" = COL_CONTROL_MOUSE, "Adn-PPARG-/-" = COL_KO_MOUSE)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    labs(title = paste0(gene_symbol_upper_safe, " expression: Adn-Pparg-KO vs Control"),
         subtitle = paste0("n = ", sum(gene_stats$n), " (", gene_stats$n[1], " per group); ", pval_label,
                            "; black points = mean \u00B1 SEM"),
         x = NULL, y = paste0("Normalized intensity (", gene_symbol_upper_safe, ")"))

  ggsave(file.path(out_dir, paste0(gene_symbol_upper_safe, "_violinplot.pdf")), p, width = 5.5, height = 6.5)
  cat("Violin plot written:", paste0(gene_symbol_upper_safe, "_violinplot.pdf"), " (", pval_label, ")\n")
}

cat("\n=== Generating individual violin plots for every single-probe panel gene ===\n")
invisible(lapply(panel_genes, make_mouse_gene_violin))
cat("(Pparg skipped above by design -- see Section 15b network for its 3 probes shown individually.)\n\n")

## ---- 13. Heatmap of the adipocyte program (Pparg excluded) ------------------

panel_probes <- manifest %>% filter(Symbol %in% panel_genes) %>% pull(Probe_Id)
panel_probes <- intersect(panel_probes, rownames(expr))
heat_mat <- log2(expr[panel_probes, ] + 1)
rownames(heat_mat) <- manifest[rownames(heat_mat), "Symbol"]
heat_z <- t(scale(t(heat_mat)))

col_anno <- data.frame(Genotype = meta$condition); rownames(col_anno) <- meta$sample
anno_colors <- list(Genotype = c("Con Adn-PPARG-/-" = "#4C72B0", "Adn-PPARG-/-" = "#DD8452"))

pdf(file.path(out_dir, "adipocyte_program_heatmap.pdf"), width = 7, height = 8)
pheatmap(heat_z, cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = col_anno, annotation_colors = anno_colors,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Adipocyte transcriptional program (Z-score)\nAdn-Pparg-KO vs Control",
         fontsize_row = 9, border_color = NA)
dev.off()
cat("Heatmap written: adipocyte_program_heatmap.pdf\n\n")

## ---- 14. Composite adipocyte-program score -----------------------------------
## Mean of standardised expression across the panel, EXCLUDING Pnpla3 itself
## (so the score isn't partly built from the gene it's being compared against).

score_genes <- setdiff(panel_genes, "Pnpla3")
score_probes <- manifest %>% filter(Symbol %in% score_genes) %>% pull(Probe_Id)
score_probes <- intersect(score_probes, rownames(expr))
score_mat <- log2(expr[score_probes, ] + 1)
score_z <- t(scale(t(score_mat)))
program_score <- colMeans(score_z, na.rm = TRUE)

score_df <- data.frame(sample = names(program_score), program_score = program_score) %>%
  left_join(meta, by = "sample") %>%
  left_join(pnpla3_vals %>% select(sample, pnpla3_expr = expr), by = "sample")

write.csv(score_df, file.path(out_dir, "adipocyte_program_score.csv"), row.names = FALSE)

p_score_group <- ggplot(score_df, aes(x = condition, y = program_score)) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(title = "Adipocyte PPARg-program score", x = NULL, y = "Mean Z-score (panel, Pnpla3 excluded)")
ggsave(file.path(out_dir, "program_score_by_group.pdf"), p_score_group, width = 5, height = 5)

p_score_vs_pnpla3 <- ggplot(score_df, aes(x = program_score, y = pnpla3_expr, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(x = program_score, y = pnpla3_expr),
              se = TRUE, color = "black", linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(title = "Pnpla3 expression vs adipocyte-program score",
       subtitle = "Interpretation: Pnpla3 covaries with the program -- not a claim that any one gene regulates it",
       x = "Adipocyte program score (mean Z, panel excl. Pnpla3)", y = "Pnpla3 (normalized intensity)", color = NULL)
ggsave(file.path(out_dir, "Pnpla3_vs_program_score.pdf"), p_score_vs_pnpla3, width = 7, height = 5.5)
cat("Composite score plots written: program_score_by_group.pdf, Pnpla3_vs_program_score.pdf\n\n")

## ---- 15. Pairwise Pnpla3 correlation scatterplots (NOT a network) -----------
## With n=6, a correlation network would be statistically fragile and easily
## misread as causal (see header). These are shown as individual scatterplots
## with all 6 points visible, framed as covariation, not regulation.

partner_genes <- c("Adipoq","Fabp4","Scd1","Cd36","Slc2a4")
partner_probes <- manifest %>% filter(Symbol %in% partner_genes) %>%
  distinct(Symbol, .keep_all = TRUE) %>%  # one probe per gene for this illustrative plot
  pull(Probe_Id)
names(partner_probes) <- manifest[partner_probes, "Symbol"]

scatter_df <- lapply(names(partner_probes), function(g) {
  data.frame(sample = colnames(expr),
             partner_gene = g,
             partner_expr = as.numeric(expr[partner_probes[g], ]),
             pnpla3_expr  = as.numeric(expr[pnpla3_probe, ]))
}) %>% bind_rows() %>% left_join(meta, by = "sample")

cor_stats <- scatter_df %>% group_by(partner_gene) %>%
  summarise(r = cor(pnpla3_expr, partner_expr, method = "pearson"),
            p = cor.test(pnpla3_expr, partner_expr)$p.value)
cat("=== Pnpla3 vs individual adipocyte genes (n=6; interpret as covariation) ===\n")
print(cor_stats)
write.csv(cor_stats, file.path(out_dir, "pnpla3_pairwise_correlations.csv"), row.names = FALSE)

p_scatter <- ggplot(scatter_df, aes(x = partner_expr, y = pnpla3_expr, color = condition)) +
  geom_point(size = 2.5) +
  facet_wrap(~partner_gene, scales = "free_x") +
  theme_bw(base_size = 10) +
  labs(title = "Pnpla3 vs individual adipocyte genes across all 6 samples",
       subtitle = "n=6 -- shown as covariation with the adipocyte program, not evidence of direct regulation",
       x = "Partner gene expression", y = "Pnpla3 expression", color = NULL)
ggsave(file.path(out_dir, "Pnpla3_pairwise_scatterplots.pdf"), p_scatter, width = 9, height = 6)
cat("Pairwise scatterplots written: Pnpla3_pairwise_scatterplots.pdf\n\n")

## ---- 15b. Gene-to-gene correlation NETWORK (node/edge graph) ----------------
## A real network graph, distinct from the pairwise scatterplots above.
## Pparg is included as 3 SEPARATE, individually labelled nodes (its 3
## Illumina probes, verified against the GSE62937 series matrix + GPL6887
## manifest: ILMN_1221060=3' end, ILMN_1221594=mid, ILMN_1221893=5' end -- see
## project notes). These probes disagree in direction and are NEVER averaged
## or summed into one "Pparg" node -- each is its own node, so nothing here
## misrepresents the underlying disagreement.

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
suppressPackageStartupMessages(library(igraph))

network_genes <- c(panel_genes, "Pparg")  # panel_genes excludes Pparg already (Section 10)
network_probes_df <- manifest %>% filter(Symbol %in% network_genes)

# Keep every probe as its own row (never collapse multi-probe genes)
network_expr <- log2(expr[network_probes_df$Probe_Id, , drop = FALSE] + 1)
node_labels <- case_when(
  network_probes_df$Symbol == "Pparg" & network_probes_df$Probe_Id == "ILMN_1221060" ~ "Pparg_3prime_probe",
  network_probes_df$Symbol == "Pparg" & network_probes_df$Probe_Id == "ILMN_1221594" ~ "Pparg_mid_probe",
  network_probes_df$Symbol == "Pparg" & network_probes_df$Probe_Id == "ILMN_1221893" ~ "Pparg_5prime_probe",
  TRUE ~ network_probes_df$Symbol
)
# Any OTHER gene with multiple probes (e.g. Lpl, Pck1 each have 2 in this panel)
# also needs a unique label -- append the probe ID rather than leaving duplicates,
# since a network graph cannot have two nodes with the same name.
dup_mask <- duplicated(node_labels) | duplicated(node_labels, fromLast = TRUE)
node_labels[dup_mask] <- paste0(node_labels[dup_mask], " (", network_probes_df$Probe_Id[dup_mask], ")")
rownames(network_expr) <- node_labels

net_cor <- cor(t(network_expr), method = "pearson")
n_net <- ncol(t(network_expr))
net_pmat <- matrix(0, nrow(net_cor), nrow(net_cor), dimnames = dimnames(net_cor))
for (i in 1:nrow(net_cor)) for (j in 1:nrow(net_cor)) {
  if (i != j) net_pmat[i, j] <- cor.test(network_expr[i, ], network_expr[j, ])$p.value
}

net_edges <- data.frame()
for (i in 1:(nrow(net_cor)-1)) for (j in (i+1):nrow(net_cor)) {
  r <- net_cor[i, j]; p <- net_pmat[i, j]
  if (abs(r) > 0.6 && p < 0.05) {
    net_edges <- rbind(net_edges, data.frame(gene1 = rownames(net_cor)[i],
                                              gene2 = rownames(net_cor)[j], r = r, p = p))
  }
}
write.csv(net_edges, file.path(out_dir, "network_significant_edges.csv"), row.names = FALSE)
cat(nrow(net_edges), "significant network edges (|r| > 0.6, p < 0.05)\n")

if (nrow(net_edges) > 0) {
  g <- graph_from_data_frame(net_edges, directed = FALSE, vertices = node_labels)
  E(g)$width <- abs(net_edges$r) * 3
  E(g)$color <- ifelse(net_edges$r > 0, COL_KO_MOUSE, COL_CONTROL_MOUSE)  # gold=+corr, indigo=-corr
  is_pparg_node <- grepl("^Pparg_", V(g)$name)
  V(g)$color <- ifelse(V(g)$name == "Pnpla3", "#D4AF37",             # muted gold for Pnpla3
                 ifelse(is_pparg_node, "#B565A7", "#CFCFCF"))         # orchid purple for Pparg probes, grey otherwise
  V(g)$size  <- ifelse(V(g)$name == "Pnpla3", 12, ifelse(is_pparg_node, 10, 8))  # smaller still -- the core cluster is a near-complete subgraph, so force-directed layouts will always pull it tight; small nodes + a large scale factor is the reliable fix

  # The densely-interconnected cluster (Pnpla3 + Pparg probes + adipogenic
  # genes) has many mutual edges, which pulls plain Fruchterman-Reingold into
  # a tight clump. Increase repulsion/iterations and then scale the resulting
  # coordinates outward to spread nodes apart while preserving the same
  # relative structure; labels are placed just below each node rather than
  # on top of it, so they stay legible once nodes are smaller and spread out.
  set.seed(42)
  lay <- layout_with_fr(g, niter = 5000,
                         area = vcount(g)^3.2,
                         repulserad = vcount(g)^3.8)
  lay <- lay * 3.5  # spread the whole layout outward substantially

  pdf(file.path(out_dir, "gene_network_graph.pdf"), width = 12, height = 12)
  plot(g, layout = lay,
       vertex.label.cex = 0.8, vertex.label.color = "black",
       vertex.label.dist = 1.5, vertex.label.degree = -pi/2,
       vertex.frame.color = "white",
       main = "Gene-gene co-expression network (|r|>0.6, p<0.05)\nGold = Pnpla3; Purple = Pparg (3 probes, shown separately)")
  dev.off()
  cat("Network graph written: gene_network_graph.pdf\n\n")
} else {
  cat("No edges passed threshold -- network not plotted.\n\n")
}

## ---- 16. Ranked log2FC plot: Pnpla3 relative to the panel -------------------

p_ranked <- ggplot(df_forest, aes(x = reorder(row_label, logFC), y = logFC)) +
  geom_col(aes(fill = Symbol == "Pnpla3")) +
  scale_fill_manual(values = c("TRUE" = "gold3", "FALSE" = "grey60"), guide = "none") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "Magnitude of Pnpla3's transcriptional dependence on Pparg,\nrelative to established adipocyte genes",
       x = NULL, y = "log2 fold change (KO vs Control)")
ggsave(file.path(out_dir, "ranked_log2FC_panel.pdf"), p_ranked, width = 7, height = 6)
cat("Ranked log2FC plot written: ranked_log2FC_panel.pdf\n\n")

## ---- 17. Observational tiers (descriptive grouping, not a statistical test) -
## Purely descriptive: groups panel genes by magnitude/significance of change.
## This is a narrative aid, not an additional statistical claim.

tiers <- panel_de %>%
  mutate(tier = case_when(
    adj.P.Val < 0.05 & abs(logFC) > 1   ~ "Strongly altered",
    adj.P.Val < 0.05 & abs(logFC) <= 1  ~ "Moderately altered",
    TRUE                                 ~ "Relatively preserved"
  )) %>%
  select(Symbol, logFC, adj.P.Val, tier) %>% arrange(tier, logFC)

write.csv(tiers, file.path(out_dir, "descriptive_tiers.csv"), row.names = FALSE)
cat("=== Descriptive tiers (narrative aid only -- not a separate statistical test) ===\n")
print(tiers)

## ---- 18. Summary ---------------------------------------------------------------

cat("\n\n===== DONE =====\n")
cat("Outputs written to:", out_dir, "\n")
if (!RUN_LIMMA) {
  cat("\n*** REMINDER: limma was not available in this session. Sections 5, 6, 7, 11\n")
  cat("*** used a Welch's-t-test fallback for demonstration. Re-run this entire\n")
  cat("*** script on a machine with limma installed before using these numbers\n")
  cat("*** in the manuscript. ***\n")
}
if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("msigdbr", quietly = TRUE)) {
  cat("\n*** REMINDER: GSEA (Section 9) requires fgsea + msigdbr -- not run. ***\n")
}
