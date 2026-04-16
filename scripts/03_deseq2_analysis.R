# 03_deseq2_analysis.R
#
# Purpose:
# Perform differential expression analysis between tumor and normal samples
# in the TCGA-COAD bulk RNA-seq dataset using DESeq2.
#
# This script:
# 1. loads the filtered SummarizedExperiment object and metadata
# 2. extracts raw count data appropriate for DESeq2
# 3. builds a DESeq2 dataset
# 4. filters low-count genes
# 5. runs differential expression analysis
# 6. saves full and significant result tables
#
# Expected outputs:
# - results/tcga_coad_deseq2_results_full.csv
# - results/tcga_coad_deseq2_results_full.rds
# - results/tcga_coad_deseq2_results_significant.csv
# - results/tcga_coad_dds.rds

# ---------------------------
# 1. Install and load packages
# ---------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_cran_packages <- c("dplyr")
for (pkg in required_cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

required_bioc_packages <- c("SummarizedExperiment", "DESeq2")
for (pkg in required_bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(SummarizedExperiment)
library(DESeq2)
library(dplyr)

# ---------------------------
# 2. Define input and output paths
# ---------------------------

input_se_path <- "data/tcga_coad_tumor_vs_normal_se.rds"
input_metadata_path <- "data/tcga_coad_metadata_tumor_vs_normal.rds"

results_full_csv_path <- "results/tcga_coad_deseq2_results_full.csv"
results_full_rds_path <- "results/tcga_coad_deseq2_results_full.rds"
results_sig_csv_path <- "results/tcga_coad_deseq2_results_significant.csv"
dds_rds_path <- "results/tcga_coad_dds.rds"

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

if (!file.exists(input_se_path)) {
  stop("Input file not found: ", input_se_path,
       "\nRun scripts/02_prepare_metadata.R first.")
}

if (!file.exists(input_metadata_path)) {
  stop("Input file not found: ", input_metadata_path,
       "\nRun scripts/02_prepare_metadata.R first.")
}

# ---------------------------
# 3. Load filtered data
# ---------------------------

coad_se_filtered <- readRDS(input_se_path)
metadata_clean <- readRDS(input_metadata_path)

# ---------------------------
# 4. Extract count matrix
# ---------------------------
# For TCGA STAR - Counts data, the SummarizedExperiment often contains
# multiple assays. DESeq2 requires raw integer counts, so we preferentially
# use the "unstranded" assay when available.

available_assays <- names(assays(coad_se_filtered))
cat("Available assays:\n")
print(available_assays)

if ("unstranded" %in% available_assays) {
  count_matrix <- assays(coad_se_filtered)[["unstranded"]]
  selected_assay <- "unstranded"
} else {
  count_matrix <- assays(coad_se_filtered)[[1]]
  selected_assay <- available_assays[1]
}

cat("\nSelected assay for DESeq2: ", selected_assay, "\n", sep = "")

# Convert to matrix explicitly
count_matrix <- as.matrix(count_matrix)

# Ensure metadata order matches count matrix columns
metadata_clean <- metadata_clean[colnames(count_matrix), , drop = FALSE]

# ---------------------------
# 5. Build DESeq2 dataset
# ---------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata_clean,
  design = ~ analysis_group
)

# ---------------------------
# 6. Filter low-count genes
# ---------------------------
# Keep genes with at least 10 counts in at least 10 samples

keep_genes <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep_genes, ]

cat("Number of genes retained after filtering: ", nrow(dds), "\n", sep = "")

# ---------------------------
# 7. Run DESeq2
# ---------------------------

dds <- DESeq(dds)

# Tumor vs Normal contrast
res <- results(
  dds,
  contrast = c("analysis_group", "Tumor", "Normal")
)

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Reorder columns
res_df <- res_df %>%
  select(gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# Order by adjusted p-value, then p-value
res_df <- res_df %>%
  arrange(padj, pvalue)

# Significant genes
res_sig_df <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05)

# ---------------------------
# 8. Save outputs
# ---------------------------

write.csv(res_df, file = results_full_csv_path, row.names = FALSE)
saveRDS(res_df, file = results_full_rds_path)
write.csv(res_sig_df, file = results_sig_csv_path, row.names = FALSE)
saveRDS(dds, file = dds_rds_path)

# ---------------------------
# 9. Print summary
# ---------------------------

cat("DESeq2 analysis complete.\n")
cat("Selected assay: ", selected_assay, "\n", sep = "")
cat("Samples analyzed: ", ncol(dds), "\n", sep = "")
cat("Genes tested after filtering: ", nrow(dds), "\n", sep = "")
cat("Significant genes (padj < 0.05): ", nrow(res_sig_df), "\n\n", sep = "")

cat("Saved files:\n")
cat("- ", results_full_csv_path, "\n", sep = "")
cat("- ", results_full_rds_path, "\n", sep = "")
cat("- ", results_sig_csv_path, "\n", sep = "")
cat("- ", dds_rds_path, "\n", sep = "")