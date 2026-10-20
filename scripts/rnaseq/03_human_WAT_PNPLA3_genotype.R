# --- Portable paths: set these two lines (or run from the project root with data/ and results/ folders) ---
DATA_DIR <- Sys.getenv("PPARG_DATA", "data"); RESULTS_DIR <- Sys.getenv("PPARG_RESULTS", "results")

#!/usr/bin/env Rscript
# ==============================================================================
# Human white adipose tissue, PNPLA3 CC vs GG (I148M): companion R pipeline
# Covers what could NOT be run in the sandbox that produced the Python-based
# results (DESeq2, GSEA Hallmark/GO, correlation network, heatmap, volcano,
# PCA) -- those are REAL, already computed and verified against your actual
# data. THIS script adds three things only R/Bioconductor can do:
#
#   1. pathview  -- draws your actual log2 fold-changes onto real KEGG pathway
#                   diagram images (e.g. "Fatty acid metabolism", "PPAR
#                   signaling pathway"). No Python equivalent was available.
#   2. clusterProfiler ORA -- the field-standard over-representation tool,
#                   complementary to the GSEA already run; reviewers will
#                   expect to see this named specifically.
#   3. WGCNA     -- with n=22 (11 vs 11), true weighted gene co-expression
#                   network module detection is now statistically defensible,
#                   unlike the n=6 mouse cohort where it was correctly avoided.
#
# IMPORTANT: this script expects the SAME filtered count matrix and metadata
# used for the Python analysis. Re-derive them exactly as follows so results
# match precisely (do not re-filter differently, or your KEGG/GO results will
# not correspond to the same gene set already reported to you).
# ==============================================================================

## ---- 0. Setup ---------------------------------------------------------------

required_bioc <- c("DESeq2", "fgsea", "clusterProfiler", "org.Hs.eg.db",
                    "pathview", "WGCNA", "AnnotationDbi")
required_cran <- c("dplyr", "tidyr", "tibble", "ggplot2", "pheatmap", "msigdbr")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(required_bioc, rownames(installed.packages()))
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
missing_cran <- setdiff(required_cran, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran)

suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(clusterProfiler); library(org.Hs.eg.db)
  library(pathview); library(WGCNA); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(pheatmap); library(msigdbr)
})

## ---- 1. Paths (EDIT IF NEEDED) ----------------------------------------------

count_path <- file.path(DATA_DIR, "countmatrix.csv")
meta_path  <- file.path(DATA_DIR, "metadata.csv")
out_dir    <- file.path(RESULTS_DIR, "human_WAT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 2. Load, aggregate duplicate symbols (SUM -- valid for RNA-seq counts,
##          unlike the earlier microarray probe datasets), filter -----------

raw <- read.csv(count_path, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)
colnames(meta) <- c("sample", "condition")
meta$condition <- factor(meta$condition, levels = c("PNPLA3 CC", "PNPLA3 GG"))
rownames(meta) <- meta$sample

sample_cols <- setdiff(colnames(raw), "gene")
gene_level <- raw %>% group_by(gene) %>% summarise(across(all_of(sample_cols), sum), .groups = "drop")

counts <- gene_level %>% column_to_rownames("gene")
counts <- round(as.matrix(counts[, meta$sample]))
mode(counts) <- "integer"

min_group_n <- min(table(meta$condition))
keep <- rowSums(counts >= 10) >= min_group_n
counts <- counts[keep, ]
cat("Genes retained after filtering:", nrow(counts), "\n")
cat("NOTE: confirm this count matches the Python analysis (11,849 genes) --\n")
cat("      if it doesn't, the filtering logic has diverged and downstream\n")
cat("      results will not correspond to what was already reported to you.\n\n")

if (!"PNPLA3" %in% rownames(counts)) {
  cat("*** PNPLA3 was filtered out (too low/sparse) -- matches the Python finding.\n")
  cat("*** This is expected and was investigated in detail already: one sample\n")
  cat("*** (P574) has an extreme outlier value; robust tests found no significant\n")
  cat("*** CC-vs-GG difference in PNPLA3 mRNA itself. See prior analysis notes.\n\n")
}

## ---- 3. DESeq2 (should match the pydeseq2 results already reported) --------

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "PNPLA3 GG", "PNPLA3 CC"))
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% arrange(padj)
write.csv(res_df, file.path(out_dir, "DESeq2_results_R.csv"), row.names = FALSE)

cat("N significant at padj<0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("(Python/pydeseq2 found 0 -- confirm this R run agrees before proceeding;\n")
cat(" small numerical differences between pydeseq2 and DESeq2 are possible but\n")
cat(" the overall conclusion -- no genes at genome-wide FDR<0.05 -- should hold.)\n\n")

## ---- 4. clusterProfiler ORA (complementary to the GSEA already run) --------
## Standard over-representation analysis on nominally significant genes
## (p<0.01, since nothing survives FDR correction here -- state this
## threshold explicitly in the manuscript as exploratory/hypothesis-generating)

sig_genes <- res_df %>% filter(pvalue < 0.01) %>% pull(gene)
cat(length(sig_genes), "genes at nominal p<0.01 (exploratory ORA input -- not FDR-significant)\n")

ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
write.csv(as.data.frame(ego), file.path(out_dir, "clusterProfiler_ORA_GO_BP.csv"), row.names = FALSE)
cat("clusterProfiler GO:BP ORA written. Top terms:\n")
print(head(as.data.frame(ego)[, c("Description", "pvalue", "p.adjust", "Count")], 15))

ekegg <- enrichKEGG(gene = bitr(sig_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID,
                     organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ekegg), file.path(out_dir, "clusterProfiler_ORA_KEGG.csv"), row.names = FALSE)
cat("\nclusterProfiler KEGG ORA written. Top pathways:\n")
print(head(as.data.frame(ekegg)[, c("Description", "pvalue", "p.adjust", "Count")], 15))

## ---- 4b. clusterProfiler plots (MISSING from the original script -- added) --
## dotplot/barplot are the standard clusterProfiler visualizations; the
## original script only wrote CSVs and console output, which was an omission.

if (nrow(as.data.frame(ego)) > 0) {
  pdf(file.path(out_dir, "clusterProfiler_GO_dotplot.pdf"), width = 9, height = 8)
  print(dotplot(ego, showCategory = 20, title = "GO Biological Process ORA (nominal p<0.01 genes)"))
  dev.off()
  cat("GO dotplot written: clusterProfiler_GO_dotplot.pdf\n")
} else {
  cat("GO ORA returned 0 enriched terms -- no dotplot produced (this itself is informative; report it).\n")
}

if (nrow(as.data.frame(ekegg)) > 0) {
  pdf(file.path(out_dir, "clusterProfiler_KEGG_barplot.pdf"), width = 9, height = 6)
  print(barplot(ekegg, showCategory = 15, title = "KEGG pathway ORA (nominal p<0.01 genes)"))
  dev.off()
  cat("KEGG barplot written: clusterProfiler_KEGG_barplot.pdf\n\n")
} else {
  cat("KEGG ORA returned 0 enriched terms -- no barplot produced (this itself is informative; report it).\n\n")
}

## ---- 5. pathview: draw actual KEGG pathway diagrams -------------------------
## This is the one analysis that genuinely cannot be reproduced outside R.
## Colors each gene on the real KEGG pathway map by its log2 fold change.
## Edit the pathway IDs below to whatever clusterProfiler's KEGG ORA (above)
## flags as significant -- these are illustrative starting points based on
## the biology already found (fatty acid metabolism, PPAR signaling,
## adipocytokine signaling, oxidative phosphorylation).

gene_data <- res_df %>% filter(!is.na(log2FoldChange)) %>%
  select(gene, log2FoldChange) %>% deframe()

entrez_map <- bitr(names(gene_data), "SYMBOL", "ENTREZID", org.Hs.eg.db)
gene_data_entrez <- gene_data[entrez_map$SYMBOL]
names(gene_data_entrez) <- entrez_map$ENTREZID

pathways_to_draw <- c(
  "hsa03320",  # PPAR signaling pathway
  "hsa00071",  # Fatty acid degradation
  "hsa04920",  # Adipocytokine signaling pathway
  "hsa00190",  # Oxidative phosphorylation
  "hsa04668"   # TNF signaling pathway
)

setwd(out_dir)  # pathview writes files to the working directory
for (pw in pathways_to_draw) {
  tryCatch({
    pathview(gene.data = gene_data_entrez, pathway.id = pw, species = "hsa",
             limit = list(gene = 2), low = "#264653", mid = "white", high = "#E76F51")
    cat("pathview diagram written for", pw, "\n")
  }, error = function(e) cat("pathview failed for", pw, ":", conditionMessage(e), "\n"))
}

## ---- 6. WGCNA: real co-expression module detection (n=22 justifies this) ---
## Unlike the n=6 mouse cohort (where WGCNA was explicitly avoided as
## statistically fragile), n=22 (11 vs 11) is within a defensible range for
## basic module detection, though still on the smaller side for WGCNA by
## field convention (ideally 20-30+ *independent* samples, which this is).

vsd <- vst(dds, blind = FALSE)
expr_mat <- t(assay(vsd))  # WGCNA wants samples x genes

# Filter to most variable genes for tractability
var_genes <- order(apply(expr_mat, 2, var), decreasing = TRUE)[1:min(3000, ncol(expr_mat))]
expr_wgcna <- expr_mat[, var_genes]

# Check for outlier samples/genes (WGCNA's own recommended first step)
gsg <- goodSamplesGenes(expr_wgcna, verbose = 3)
if (!gsg$allOK) expr_wgcna <- expr_wgcna[gsg$goodSamples, gsg$goodGenes]
cat("WGCNA input matrix after goodSamplesGenes check:", nrow(expr_wgcna), "samples x",
    ncol(expr_wgcna), "genes\n\n")

powers <- c(1:20)
sft <- tryCatch(
  pickSoftThreshold(expr_wgcna, powerVector = powers, verbose = 0),
  error = function(e) { cat("*** pickSoftThreshold FAILED:", conditionMessage(e), "***\n"); NULL }
)

if (!is.null(sft)) {
  tryCatch({
    pdf(file.path(out_dir, "WGCNA_soft_threshold.pdf"), width = 8, height = 5)
    par(mfrow = c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)", ylab="Scale Free Topology Fit R^2", main="Scale independence")
    abline(h=0.85, col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)", ylab="Mean Connectivity", main="Mean connectivity")
    dev.off()
    cat("WGCNA_soft_threshold.pdf written successfully.\n")
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()  # ensure the graphics device is never left open on failure
    cat("*** Soft-threshold PLOT failed:", conditionMessage(e), "***\n")
  })
  cat("\nCHOOSE the power where the R^2 curve first plateaus above ~0.85 (read off the plot)\n")
  cat("and set SOFT_POWER below before trusting downstream modules; do not use a default blindly.\n\n")
} else {
  cat("*** Skipping soft-threshold plot -- pickSoftThreshold did not return a result. ***\n")
  cat("*** This can happen with small/borderline sample sizes; check the console output ***\n")
  cat("*** above this point for WGCNA's own warnings before proceeding further.        ***\n\n")
}

SOFT_POWER <- 6  # <-- EDIT based on the plot above before trusting downstream modules

net <- tryCatch(
  blockwiseModules(expr_wgcna, power = SOFT_POWER, TOMType = "signed",
                    minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25,
                    numericLabels = FALSE, saveTOMs = FALSE, verbose = 3),
  error = function(e) { cat("*** blockwiseModules FAILED:", conditionMessage(e), "***\n"); NULL }
)

if (!is.null(net)) {
  module_membership <- data.frame(gene = colnames(expr_wgcna), module = net$colors)
  write.csv(module_membership, file.path(out_dir, "WGCNA_module_membership.csv"), row.names = FALSE)
  cat("WGCNA_module_membership.csv written.\n")

  cat("\nModule sizes:\n")
  print(table(net$colors))

  if ("PNPLA3" %in% module_membership$gene) {
    cat("\nPNPLA3's module:", module_membership$module[module_membership$gene == "PNPLA3"], "\n")
  } else {
    cat("\nPNPLA3 not among the most-variable genes used for WGCNA (consistent with\n")
    cat("its low/sparse expression noted earlier) -- module membership not available for it.\n")
  }
  if ("PPARG" %in% module_membership$gene) {
    cat("PPARG's module:", module_membership$module[module_membership$gene == "PPARG"], "\n")
  }

  tryCatch({
    pdf(file.path(out_dir, "WGCNA_dendrogram.pdf"), width = 10, height = 6)
    plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                         "Module colors", dendroLabels = FALSE, hang = 0.03,
                         addGuide = TRUE, guideHang = 0.05)
    dev.off()
    cat("WGCNA_dendrogram.pdf written successfully.\n")
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    cat("*** Dendrogram PLOT failed:", conditionMessage(e), "***\n")
  })
} else {
  cat("*** Skipping module detection and dendrogram -- blockwiseModules did not complete. ***\n")
  cat("*** Common causes at this sample size: too few samples for the chosen power, or   ***\n")
  cat("*** insufficient memory for 3000 genes -- try reducing to var_genes <- ...[1:1000] ***\n")
  cat("*** near the top of this section and re-running just Section 6.                   ***\n")
}

cat("\n===== DONE =====\n")
cat("Outputs written to:", out_dir, "\n")
cat(" - DESeq2_results_R.csv (should match the Python pydeseq2 results already reported)\n")
cat(" - clusterProfiler_ORA_GO_BP.csv / clusterProfiler_ORA_KEGG.csv\n")
cat(" - KEGG pathway diagrams (hsa*.pathview.png files, in out_dir)\n")
cat(" - WGCNA_soft_threshold.pdf, WGCNA_dendrogram.pdf, WGCNA_module_membership.csv\n")

## ==============================================================================
## ---- 7. GSVA/ssGSEA per-sample pathway scores (already run and verified in
##          Python -- this R version lets you cite the field-standard GSVA
##          package by name and reproduce it locally)
## ==============================================================================

if (!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA", update = FALSE, ask = FALSE)
library(GSVA)

vsd_mat <- assay(vst(dds, blind = FALSE))

hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
target_sets <- hallmark_sets[c("HALLMARK_FATTY_ACID_METABOLISM",
                                "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                "HALLMARK_TNFA_SIGNALING_VIA_NFKB")]

gsvaPar <- ssgseaParam(vsd_mat, target_sets)
gsva_scores <- gsva(gsvaPar)
gsva_df <- as.data.frame(t(gsva_scores)) %>% rownames_to_column("sample") %>%
  left_join(meta, by = "sample")
write.csv(gsva_df, file.path(out_dir, "GSVA_ssGSEA_scores_R.csv"), row.names = FALSE)
cat("\nGSVA/ssGSEA scores written: GSVA_ssGSEA_scores_R.csv\n")
cat("(Compare against the Python ssGSEA_scores_human_WAT.csv already reported --\n")
cat(" values will differ slightly by implementation but direction should agree.)\n\n")

## ==============================================================================
## ---- 8. Targeted co-expression network (Hmisc::rcorr) -- NOT WGCNA
## As recommended: a focused, hypothesis-driven co-expression network on the
## specific FA-oxidation and inflammatory gene panels, rather than a
## genome-wide or WGCNA-based network. This is a co-expression network, not a
## regulatory network -- do not describe it as showing what "regulates" what.
## ==============================================================================

if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(Hmisc); library(igraph)

fa_genes <- c("ACOX1","CPT1A","CPT2","ACADVL","HADHA","HADHB","ECHS1","ETFDH",
              "NDUFS1","NDUFV1","SDHA","SDHB","UQCRC1","COX4I1")
inflam_genes <- c("CCL2","NFKB2","CXCL2","DUSP5","HES1","MRC1","CD68","ITGAM","TLR4")
targeted_genes <- intersect(c(fa_genes, inflam_genes), rownames(vsd_mat))

targeted_mat <- t(vsd_mat[targeted_genes, ])
rcorr_res <- rcorr(targeted_mat, type = "spearman")
r_mat <- rcorr_res$r
p_mat <- rcorr_res$P
p_mat[is.na(p_mat)] <- 1
p_adj <- matrix(p.adjust(p_mat, method = "BH"), nrow = nrow(p_mat), dimnames = dimnames(p_mat))

edges <- data.frame()
for (i in 1:(length(targeted_genes)-1)) for (j in (i+1):length(targeted_genes)) {
  r <- r_mat[i,j]; padj <- p_adj[i,j]
  if (abs(r) > 0.5 && padj < 0.05) {
    edges <- rbind(edges, data.frame(gene1 = targeted_genes[i], gene2 = targeted_genes[j], r = r, padj = padj))
  }
}
write.csv(edges, file.path(out_dir, "targeted_coexpression_edges.csv"), row.names = FALSE)
cat(nrow(edges), "significant edges (|rho|>0.5, BH-adjusted p<0.05) in the targeted co-expression network\n")

if (nrow(edges) > 0) {
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = targeted_genes)
  V(g)$group <- ifelse(V(g)$name %in% fa_genes, "FA oxidation/OXPHOS", "Inflammatory")
  V(g)$color <- ifelse(V(g)$group == "FA oxidation/OXPHOS", "#264653", "#E76F51")
  E(g)$color <- ifelse(edges$r > 0, "grey40", "grey70")
  pdf(file.path(out_dir, "targeted_coexpression_network.pdf"), width = 9, height = 9)
  plot(g, layout = layout_with_fr(g), vertex.size = 14, vertex.label.cex = 0.8,
       vertex.label.color = "black", edge.width = abs(edges$r)*3,
       main = "Targeted co-expression network (Spearman, BH-adjusted p<0.05)\nTeal = FA oxidation/OXPHOS, Coral = Inflammatory\n(co-expression only -- not a regulatory network)")
  dev.off()
  cat("targeted_coexpression_network.pdf written.\n")
}

## ==============================================================================
## ---- 9. Transcript/isoform-level analysis -- CONDITIONAL, requires raw FASTQ
## Gene-level counts (as used throughout this analysis) cannot answer whether
## genotype affects PNPLA3 transcript/isoform USAGE -- that requires
## transcript-level quantification from raw reads. This section is a template
## only; it will not run without your own Salmon/kallisto output.
##
##   FASTQ -> Salmon/kallisto (transcript-level quantification)
##          -> tximport (import into R)
##          -> IsoformSwitchAnalyzeR or DRIMSeq/DEXSeq (differential usage)
##
## Do NOT conflate this with allele-specific expression (ASE): ASE requires
## heterozygous (CG) individuals to compare reference vs alternate allele
## read counts at rs738409 WITHIN each person. Your cohort (CC/GG homozygotes
## only) cannot support an ASE analysis -- that would need a separate cohort
## with CG heterozygotes and allele-aware read counting (e.g. WASP, phASER).
## ==============================================================================

cat("\n*** Transcript-level / isoform-usage analysis requires raw FASTQ files and \n")
cat("*** Salmon/kallisto + tximport + IsoformSwitchAnalyzeR -- not run here since \n")
cat("*** only gene-level counts were provided. See comments above for the pipeline \n")
cat("*** if you obtain transcript-level quantification. \n")
cat("*** Allele-specific expression is a DIFFERENT analysis requiring CG heterozygotes \n")
cat("*** -- not possible with a CC/GG-only cohort. \n\n")

cat("===== FULL PIPELINE DONE =====\n")
