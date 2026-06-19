# ============================================================
# PHASE 2 — DESeq2 ANALYSIS FOR TP53 STATUS
# TCGA-BRCA Transcriptomic Dysregulation
# ============================================================

# OBJECTIVE
# ----------
# Compare:
#   TP53-mutant vs TP53-wildtype breast tumors
#
# OUTPUTS
# --------
# 1. Normalized counts
# 2. PCA plot
# 3. DEG table
# 4. Significant DEG table
# 5. Volcano plot
# 6. Heatmap of top DEGs
#
# ============================================================


# ============================================================
# STEP 1 — LOAD LIBRARIES
# ============================================================

required_packages <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "EnhancedVolcano",
  "dplyr",
  "RColorBrewer"
)

for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if(!requireNamespace("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
}

if(!requireNamespace("EnhancedVolcano", quietly = TRUE)){
  BiocManager::install("EnhancedVolcano")
}


# ============================================================
# STEP 2 — DEFINE PATHS
# ============================================================

base_dir <- "D:/TP53_project/data/expression/"

counts_path <- paste0(base_dir, "deseq2_ready_counts.csv.gz")
metadata_path <- paste0(base_dir, "deseq2_ready_metadata.csv")

output_dir <- paste0(base_dir, "DEA_results/")

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}


# ============================================================
# STEP 3 — LOAD DATA
# ============================================================

cat("\nLoading count matrix...\n")

counts <- read.csv(
  gzfile(counts_path),
  row.names = 1,
  check.names = FALSE
)

cat("Count matrix dimensions:\n")
print(dim(counts))

cat("\nLoading metadata...\n")

metadata <- read.csv(
  metadata_path,
  row.names = 1
)

cat("Metadata dimensions:\n")
print(dim(metadata))


# ============================================================
# STEP 4 — ENSURE MATCHING SAMPLE ORDER
# ============================================================

cat("\nChecking sample consistency...\n")

common_samples <- intersect(
  colnames(counts),
  rownames(metadata)
)

counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

cat("Matched samples:", length(common_samples), "\n")


# ============================================================
# STEP 5 — CLEAN COUNT MATRIX
# ============================================================

cat("\nCleaning count matrix...\n")

# Remove duplicated genes
counts <- counts[!duplicated(rownames(counts)), ]

# Remove low-count genes
keep <- rowSums(counts >= 10) >= 5
counts <- counts[keep, ]

cat("Remaining genes after filtering:\n")
print(nrow(counts))


# ============================================================
# STEP 6 — FACTOR SETUP
# ============================================================

metadata$TP53_status <- factor(
  metadata$TP53_status,
  levels = c("TP53_wildtype", "TP53_mutant")
)

table(metadata$TP53_status)


# ============================================================
# STEP 7 — CREATE DESEQ2 OBJECT
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = metadata,
  design = ~ TP53_status
)

cat("\nDESeq2 object created.\n")


# ============================================================
# STEP 8 — RUN DESEQ2
# ============================================================

cat("\nRunning DESeq2 analysis...\n")

dds <- DESeq(dds)

cat("DESeq2 completed.\n")


# ============================================================
# STEP 9 — NORMALIZED COUNTS
# ============================================================

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(
  normalized_counts,
  paste0(output_dir, "normalized_counts.csv")
)

cat("\nNormalized counts saved.\n")


# ============================================================
# STEP 10 — PCA PLOT
# ============================================================

cat("\nGenerating PCA plot...\n")

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = "TP53_status",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, color = TP53_status)
) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Figure 4: PCA Plot of TP53 Mutation Status") +
  theme_bw(base_size = 14)

ggsave(
  paste0(output_dir, "Figure4_PCA_plot.png"),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved.\n")


# ============================================================
# STEP 11 — DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

cat("\nExtracting DEGs...\n")

res <- results(dds)

res <- lfcShrink(
  dds,
  coef = "TP53_status_TP53_mutant_vs_TP53_wildtype",
  type = "apeglm"
)

res_df <- as.data.frame(res)

res_df <- res_df[order(res_df$padj), ]

write.csv(
  res_df,
  paste0(output_dir, "DESeq2_all_results.csv")
)

cat("Full DEG table saved.\n")


# ============================================================
# STEP 12 — FILTER SIGNIFICANT DEGs
# ============================================================

sig_res <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

cat("\nSignificant DEGs:\n")
print(nrow(sig_res))

write.csv(
  sig_res,
  paste0(output_dir, "DESeq2_significant_DEGs.csv")
)

cat("Significant DEG table saved.\n")


# ============================================================
# STEP 13 — UP/DOWN REGULATED GENES
# ============================================================

upregulated <- subset(
  sig_res,
  log2FoldChange > 1
)

downregulated <- subset(
  sig_res,
  log2FoldChange < -1
)

write.csv(
  upregulated,
  paste0(output_dir, "upregulated_genes.csv")
)

write.csv(
  downregulated,
  paste0(output_dir, "downregulated_genes.csv")
)

cat("\nUpregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")


# ============================================================
# STEP 14 — VOLCANO PLOT
# ============================================================

cat("\nGenerating volcano plot...\n")

png(
  filename = paste0(output_dir, "Figure5_Volcano_plot.png"),
  width = 3000,
  height = 2400,
  res = 300
)

EnhancedVolcano(
  res_df,
  lab = NA,   # removes all point labels (Ensembl IDs)
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 3.5,
  title = "Figure 5: TP53 Mutation-Associated DEGs",
  subtitle = "TCGA-BRCA",
  caption = paste0(
    "Total significant genes: ",
    nrow(sig_res)
  )
)
dev.off()

cat("Volcano plot saved.\n")


# ============================================================
# STEP 15 — HEATMAP OF TOP DEGs
# ============================================================

cat("\nGenerating heatmap...\n")

top_genes <- rownames(sig_res)[1:50]

mat <- assay(vsd)[top_genes, ]

mat <- t(scale(t(mat)))

annotation_col <- data.frame(
  TP53_status = metadata$TP53_status
)

rownames(annotation_col) <- rownames(metadata)

png(
  filename = paste0(output_dir, "Figure6_Heatmap_top_DEGs.png"),
  width = 3000,
  height = 3600,
  res = 300
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Figure 6: Top Differentially Expressed Genes"
)

dev.off()

cat("Heatmap saved.\n")


# ============================================================
# STEP 16 — SUMMARY REPORT
# ============================================================

summary_file <- paste0(output_dir, "analysis_summary.txt")

sink(summary_file)

cat("TP53 Differential Expression Analysis Summary\n")
cat("===========================================\n\n")

cat("Total samples:\n")
print(table(metadata$TP53_status))

cat("\nTotal genes analyzed:\n")
print(nrow(res_df))

cat("\nSignificant DEGs:\n")
print(nrow(sig_res))

cat("\nUpregulated genes:\n")
print(nrow(upregulated))

cat("\nDownregulated genes:\n")
print(nrow(downregulated))

sink()

cat("\nSummary report saved.\n")


# ============================================================
# FINAL MESSAGE
# ============================================================

cat("\n============================================\n")
cat("DESeq2 ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("\nResults directory:\n")
cat(output_dir)

cat("\n\nGenerated outputs:\n")
cat("- PCA plot\n")
cat("- Volcano plot\n")
cat("- Heatmap\n")
cat("- DEG tables\n")
cat("- Normalized counts\n")
cat("- Summary report\n")

cat("\n============================================\n")


# ============================================================
# STEP 15 — ENSEMBL → ENTREZ + GENE SYMBOL CONVERSION
# ============================================================

# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "dplyr",
  "readr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# ============================================================
# INPUT / OUTPUT PATHS
# ============================================================

input_dir  <- "D:/TP53_project/data/expression/DEA_results/"
output_dir <- "D:/TP53_project/data/expression/DEA_results/"

# ============================================================
# LOAD FILES
# ============================================================

cat("\nLoading DEG files...\n")

all_deg  <- read.csv(
  paste0(input_dir, "DESeq2_all_results.csv"),
  row.names = 1
)

sig_deg  <- read.csv(
  paste0(input_dir, "DESeq2_significant_DEGs.csv"),
  row.names = 1
)

up_deg <- read.csv(
  paste0(input_dir, "upregulated_genes.csv"),
  row.names = 1
)

down_deg <- read.csv(
  paste0(input_dir, "downregulated_genes.csv"),
  row.names = 1
)

cat("Files loaded successfully.\n")

# ============================================================
# FUNCTION — ENSEMBL CONVERSION
# ============================================================

convert_ensembl <- function(df) {
  
  # Extract Ensembl IDs from rownames
  df$ENSEMBL <- rownames(df)
  
  # Remove version numbers if present
  df$ENSEMBL <- gsub("\\..*", "", df$ENSEMBL)
  
  # Convert IDs
  gene_annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = df$ENSEMBL,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "ENSEMBL"
  )
  
  # Remove duplicated mappings
  gene_annotation <- gene_annotation %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  # Merge with DEG dataframe
  merged_df <- merge(
    df,
    gene_annotation,
    by = "ENSEMBL",
    all.x = TRUE
  )
  
  # Reorder columns
  important_cols <- c(
    "ENSEMBL",
    "SYMBOL",
    "ENTREZID"
  )
  
  other_cols <- setdiff(colnames(merged_df), important_cols)
  
  merged_df <- merged_df[, c(important_cols, other_cols)]
  
  return(merged_df)
}

# ============================================================
# CONVERT ALL FILES
# ============================================================

cat("\nConverting ENSEMBL IDs...\n")

all_deg_annotated  <- convert_ensembl(all_deg)
sig_deg_annotated  <- convert_ensembl(sig_deg)
up_deg_annotated   <- convert_ensembl(up_deg)
down_deg_annotated <- convert_ensembl(down_deg)

cat("Conversion complete.\n")

# ============================================================
# QC REPORT
# ============================================================

cat("\n================ QC REPORT ================\n")

cat("\nAll genes mapped:\n")
cat(sum(!is.na(all_deg_annotated$SYMBOL)), "\n")

cat("\nSignificant genes mapped:\n")
cat(sum(!is.na(sig_deg_annotated$SYMBOL)), "\n")

cat("\nUpregulated genes mapped:\n")
cat(sum(!is.na(up_deg_annotated$SYMBOL)), "\n")

cat("\nDownregulated genes mapped:\n")
cat(sum(!is.na(down_deg_annotated$SYMBOL)), "\n")

# ============================================================
# SAVE OUTPUTS
# ============================================================

write.csv(
  all_deg_annotated,
  paste0(output_dir, "DESeq2_all_results_annotated.csv"),
  row.names = FALSE
)

write.csv(
  sig_deg_annotated,
  paste0(output_dir, "DESeq2_significant_DEGs_annotated.csv"),
  row.names = FALSE
)

write.csv(
  up_deg_annotated,
  paste0(output_dir, "upregulated_genes_annotated.csv"),
  row.names = FALSE
)

write.csv(
  down_deg_annotated,
  paste0(output_dir, "downregulated_genes_annotated.csv"),
  row.names = FALSE
)

cat("\n===========================================\n")
cat("ANNOTATED FILES SAVED SUCCESSFULLY\n")
cat("===========================================\n")

cat("\nSaved files:\n")

cat("1. DESeq2_all_results_annotated.csv\n")
cat("2. DESeq2_significant_DEGs_annotated.csv\n")
cat("3. upregulated_genes_annotated.csv\n")
cat("4. downregulated_genes_annotated.csv\n")

# ============================================================
# STEP 16 — FUNCTIONAL ENRICHMENT ANALYSIS (FEA)
# TP53 Mutation-associated Transcriptomic Dysregulation
# ============================================================

# ============================================================
# REQUIRED PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2",
  "dplyr",
  "readr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  
  library(pkg, character.only = TRUE)
}

# ============================================================
# INPUT / OUTPUT PATHS
# ============================================================

input_dir  <- "D:/TP53_project/data/expression/DEA_results/"
output_dir <- "D:/TP53_project/data/expression/FEA_results/"

dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# LOAD ANNOTATED DEG FILE
# ============================================================

cat("\nLoading annotated DEG file...\n")

deg_df <- read.csv(
  paste0(input_dir, "DESeq2_significant_DEGs_annotated.csv")
)

cat("Rows loaded:", nrow(deg_df), "\n")

# ============================================================
# CLEAN DATA
# ============================================================

deg_df <- deg_df %>%
  filter(!is.na(ENTREZID))

deg_df$ENTREZID <- as.character(deg_df$ENTREZID)

cat("Genes with valid ENTREZ IDs:", nrow(deg_df), "\n")

# ============================================================
# SPLIT INTO UP / DOWN
# ============================================================

up_genes <- deg_df %>%
  filter(log2FoldChange > 1,
         padj < 0.05) %>%
  pull(ENTREZID)

down_genes <- deg_df %>%
  filter(log2FoldChange < -1,
         padj < 0.05) %>%
  pull(ENTREZID)

all_genes <- unique(deg_df$ENTREZID)

cat("\nUpregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")

# ============================================================
# GO ENRICHMENT — BIOLOGICAL PROCESS
# ============================================================

cat("\nRunning GO-BP enrichment...\n")

go_bp <- enrichGO(
  gene          = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ============================================================
# GO ENRICHMENT — MOLECULAR FUNCTION
# ============================================================

cat("Running GO-MF enrichment...\n")

go_mf <- enrichGO(
  gene          = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ============================================================
# GO ENRICHMENT — CELLULAR COMPONENT
# ============================================================

cat("Running GO-CC enrichment...\n")

go_cc <- enrichGO(
  gene          = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ============================================================
# KEGG ENRICHMENT
# ============================================================

cat("Running KEGG enrichment...\n")

kegg_res <- enrichKEGG(
  gene         = all_genes,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# ============================================================
# REACTOME ENRICHMENT
# ============================================================

cat("Running Reactome enrichment...\n")

reactome_res <- enrichPathway(
  gene          = all_genes,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

# ============================================================
# SAVE ENRICHMENT TABLES
# ============================================================

cat("\nSaving enrichment tables...\n")

write.csv(
  as.data.frame(go_bp),
  paste0(output_dir, "GO_BP_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(go_mf),
  paste0(output_dir, "GO_MF_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(go_cc),
  paste0(output_dir, "GO_CC_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(kegg_res),
  paste0(output_dir, "KEGG_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(reactome_res),
  paste0(output_dir, "Reactome_results.csv"),
  row.names = FALSE
)

# ============================================================
# GLOBAL PLOT STYLE
# ============================================================

theme_set(theme_bw(base_size = 14))

# ============================================================
# FIGURE 7A — GO BIOLOGICAL PROCESS DOTPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7A_GO_BP_dotplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  go_bp,
  showCategory = 15,
  title = "GO Biological Process Enrichment"
)

dev.off()

# ============================================================
# FIGURE 7B — KEGG DOTPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7B_KEGG_dotplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  kegg_res,
  showCategory = 15,
  title = "KEGG Pathway Enrichment"
)

dev.off()

# ============================================================
# FIGURE 7C — REACTOME DOTPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7C_Reactome_dotplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  reactome_res,
  showCategory = 15,
  title = "Reactome Pathway Enrichment"
)

dev.off()

# ============================================================
# FIGURE 7D — GO BARPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7D_GO_BP_barplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

barplot(
  go_bp,
  showCategory = 15,
  title = "GO Biological Process Enrichment"
)

dev.off()

# ============================================================
# FIGURE 7E — KEGG BARPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7E_KEGG_barplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

barplot(
  kegg_res,
  showCategory = 15,
  title = "KEGG Pathway Enrichment"
)

dev.off()

# ============================================================
# FIGURE 7F — REACTOME BARPLOT
# ============================================================

png(
  filename = paste0(output_dir, "Figure7F_Reactome_barplot.png"),
  width = 3200,
  height = 2400,
  res = 300
)

barplot(
  reactome_res,
  showCategory = 15,
  title = "Reactome Pathway Enrichment"
)

dev.off()

# ============================================================
# SUMMARY REPORT
# ============================================================

cat("\n====================================================\n")
cat("FUNCTIONAL ENRICHMENT ANALYSIS COMPLETE\n")
cat("====================================================\n")

cat("\nGO-BP pathways:", nrow(as.data.frame(go_bp)), "\n")
cat("GO-MF pathways:", nrow(as.data.frame(go_mf)), "\n")
cat("GO-CC pathways:", nrow(as.data.frame(go_cc)), "\n")
cat("KEGG pathways :", nrow(as.data.frame(kegg_res)), "\n")
cat("Reactome pathways:", nrow(as.data.frame(reactome_res)), "\n")

cat("\nResults directory:\n")
cat(output_dir, "\n")

cat("\nGenerated figures:\n")
cat("Figure7A_GO_BP_dotplot.png\n")
cat("Figure7B_KEGG_dotplot.png\n")
cat("Figure7C_Reactome_dotplot.png\n")
cat("Figure7D_GO_BP_barplot.png\n")
cat("Figure7E_KEGG_barplot.png\n")
cat("Figure7F_Reactome_barplot.png\n")

cat("\nCSV enrichment tables saved successfully.\n")

###############################
# ============================================================
# LOAD ANNOTATED SIGNIFICANT DEGs
# ============================================================

library(dplyr)

sig_deg_annotated <- read.csv(
  "D:/TP53_project/data/expression/DEA_results/DESeq2_significant_DEGs_annotated.csv"
)

cat("File loaded successfully.\n")

# ============================================================
# FILTER GENES FOR PPI NETWORK
# ============================================================

ppi_genes <- sig_deg_annotated %>%
  filter(
    padj < 0.01,
    abs(log2FoldChange) > 1.5,
    !is.na(SYMBOL)
  )

cat("\nGenes selected for PPI:\n")
print(nrow(ppi_genes))

# ============================================================
# EXPORT SYMBOL LIST FOR STRING
# ============================================================

write.table(
  unique(ppi_genes$SYMBOL),
  file = "D:/TP53_project/data/network_analysis/PPI_gene_symbols.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\nPPI gene list exported successfully.\n")