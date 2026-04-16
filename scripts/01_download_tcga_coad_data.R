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

# Convert list columns to character so they can be written to CSV
flatten_list_column <- function(x) {
  vapply(x, function(y) {
    if (is.null(y) || length(y) == 0) {
      return(NA_character_)
    }
    paste(as.character(y), collapse = "; ")
  }, character(1))
}

list_cols <- vapply(metadata_df, is.list, logical(1))
if (any(list_cols)) {
  metadata_df[list_cols] <- lapply(metadata_df[list_cols], flatten_list_column)
}

write.csv(metadata_df, file = "data/tcga_coad_metadata.csv", row.names = FALSE)
cat("Download and preparation complete.\n")
cat("Files saved:\n")
cat("- data/tcga_coad_summarized_experiment.rds\n")
cat("- data/tcga_coad_metadata.csv\n")
