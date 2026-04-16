# 01_download_tcga_coad_data.R
# Download TCGA-COAD RNA-seq count data and metadata using TCGAbiolinks

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages if needed
required_packages <- c("TCGAbiolinks", "SummarizedExperiment")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)

# Create data directory if it does not exist
if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

# Query TCGA-COAD RNA-seq count data
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download data
GDCdownload(query)

# Prepare SummarizedExperiment object
coad_se <- GDCprepare(query)

# Save the complete object
saveRDS(coad_se, file = "data/tcga_coad_summarized_experiment.rds")

# Save metadata separately for easier inspection
metadata_df <- as.data.frame(colData(coad_se))
write.csv(metadata_df, file = "data/tcga_coad_metadata.csv", row.names = FALSE)

cat("Download and preparation complete.\n")
cat("Files saved:\n")
cat("- data/tcga_coad_summarized_experiment.rds\n")
cat("- data/tcga_coad_metadata.csv\n")
