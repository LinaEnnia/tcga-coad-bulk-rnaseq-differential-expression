# 02_prepare_metadata.R
#
# Purpose:
# Prepare a clean tumor-vs-normal sample metadata table for the TCGA-COAD
# bulk RNA-seq differential expression analysis.
#
# This script:
# 1. loads the downloaded SummarizedExperiment object
# 2. derives TCGA sample type information from sample barcodes
# 3. keeps only tumor and normal samples relevant for the analysis
# 4. saves a compact metadata table and a filtered SummarizedExperiment object
#
# Expected outputs:
# - data/tcga_coad_metadata_tumor_vs_normal.csv
# - data/tcga_coad_metadata_tumor_vs_normal.rds
# - data/tcga_coad_tumor_vs_normal_se.rds

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

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment", ask = FALSE, update = FALSE)
}

library(SummarizedExperiment)
library(dplyr)

# ---------------------------
# 2. Define input and output paths
# ---------------------------

input_path <- "data/tcga_coad_summarized_experiment.rds"
metadata_csv_path <- "data/tcga_coad_metadata_tumor_vs_normal.csv"
metadata_rds_path <- "data/tcga_coad_metadata_tumor_vs_normal.rds"
filtered_se_path <- "data/tcga_coad_tumor_vs_normal_se.rds"

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path,
       "\nRun scripts/01_download_tcga_coad_data.R first.")
}

# ---------------------------
# 3. Load SummarizedExperiment
# ---------------------------

coad_se <- readRDS(input_path)

# ---------------------------
# 4. Build compact metadata table
# ---------------------------
# TCGA sample type can be derived from barcode characters 14-15:
# - 01 = Primary Solid Tumor
# - 11 = Solid Tissue Normal

sample_barcodes <- colnames(coad_se)
sample_type_code <- substr(sample_barcodes, 14, 15)

sample_type_label <- dplyr::case_when(
  sample_type_code == "01" ~ "Primary Solid Tumor",
  sample_type_code == "11" ~ "Solid Tissue Normal",
  TRUE ~ "Other"
)

analysis_group <- dplyr::case_when(
  sample_type_code == "01" ~ "Tumor",
  sample_type_code == "11" ~ "Normal",
  TRUE ~ NA_character_
)

metadata_clean <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id = substr(sample_barcodes, 1, 12),
  sample_type_code = sample_type_code,
  sample_type_label = sample_type_label,
  analysis_group = analysis_group,
  stringsAsFactors = FALSE
)

# Keep only tumor and normal samples for the DE analysis
metadata_clean <- metadata_clean[!is.na(metadata_clean$analysis_group), , drop = FALSE]

# Set row names to sample barcodes for downstream matching
rownames(metadata_clean) <- metadata_clean$sample_barcode

# Use a consistent reference order for DESeq2 later
metadata_clean$analysis_group <- factor(
  metadata_clean$analysis_group,
  levels = c("Normal", "Tumor")
)

# ---------------------------
# 5. Filter the SummarizedExperiment
# ---------------------------

coad_se_filtered <- coad_se[, metadata_clean$sample_barcode]

# ---------------------------
# 6. Save outputs
# ---------------------------

write.csv(metadata_clean, file = metadata_csv_path, row.names = FALSE)
saveRDS(metadata_clean, file = metadata_rds_path)
saveRDS(coad_se_filtered, file = filtered_se_path)

# ---------------------------
# 7. Print summary
# ---------------------------

cat("Metadata preparation complete.\n")
cat("Samples retained for tumor-vs-normal analysis:\n")
print(table(metadata_clean$analysis_group))

cat("\nSaved files:\n")
cat("- ", metadata_csv_path, "\n", sep = "")
cat("- ", metadata_rds_path, "\n", sep = "")
cat("- ", filtered_se_path, "\n", sep = "")