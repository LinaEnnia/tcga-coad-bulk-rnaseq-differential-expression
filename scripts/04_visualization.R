# 04_visualization.R
#
# Purpose:
# Generate key visualizations for the TCGA-COAD tumor-vs-normal
# differential expression analysis.
#
# This script:
# 1. loads the DESeq2 object and full results table
# 2. adds gene symbol annotation
# 3. applies variance stabilizing transformation (VST)
# 4. creates a PCA plot
# 5. creates a volcano plot
# 6. creates an MA plot
# 7. creates a heatmap of the top 30 most significant genes
#
# Expected outputs:
# - figures/tcga_coad_pca_plot.png
# - figures/tcga_coad_volcano_plot.png
# - figures/tcga_coad_ma_plot.png
# - figures/tcga_coad_top30_heatmap.png
# - results/tcga_coad_top30_heatmap_genes.csv

# ---------------------------
# 1. Install and load packages
# ---------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_cran_packages <- c("dplyr", "pheatmap")
for (pkg in required_cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

required_bioc_packages <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db")
for (pkg in required_bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(DESeq2)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ---------------------------
# 2. Define input and output paths
# ---------------------------

dds_rds_path <- "results/tcga_coad_dds.rds"
results_full_rds_path <- "results/tcga_coad_deseq2_results_full.rds"

pca_plot_path <- "figures/tcga_coad_pca_plot.png"
volcano_plot_path <- "figures/tcga_coad_volcano_plot.png"
ma_plot_path <- "figures/tcga_coad_ma_plot.png"
heatmap_plot_path <- "figures/tcga_coad_top30_heatmap.png"
top_genes_csv_path <- "results/tcga_coad_top30_heatmap_genes.csv"

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

if (!file.exists(dds_rds_path)) {
  stop("Input file not found: ", dds_rds_path,
       "\nRun scripts/03_deseq2_analysis.R first.")
}

if (!file.exists(results_full_rds_path)) {
  stop("Input file not found: ", results_full_rds_path,
       "\nRun scripts/03_deseq2_analysis.R first.")
}

# ---------------------------
# 3. Load DESeq2 object and results
# ---------------------------

dds <- readRDS(dds_rds_path)
res_df <- readRDS(results_full_rds_path)

# Recreate DESeq2 results object for MA plot
res <- results(
  dds,
  contrast = c("analysis_group", "Tumor", "Normal")
)

# ---------------------------
# 3b. Add gene symbol annotation
# ---------------------------
# Remove Ensembl version suffix (for example ENSG000001234.5 -> ENSG000001234)

res_df$ensembl_id_clean <- sub("\\..*$", "", res_df$gene_id)

gene_symbol_map <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = unique(res_df$ensembl_id_clean),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df$gene_symbol <- unname(gene_symbol_map[res_df$ensembl_id_clean])

# Use gene symbol when available; otherwise fall back to Ensembl ID
res_df$plot_label <- ifelse(
  !is.na(res_df$gene_symbol) & res_df$gene_symbol != "",
  res_df$gene_symbol,
  res_df$gene_id
)

# ---------------------------
# 4. Variance stabilizing transformation
# ---------------------------

vsd <- vst(dds, blind = FALSE)

# ---------------------------
# 5. PCA plot
# ---------------------------

pca_matrix <- assay(vsd)

pca_res <- prcomp(t(pca_matrix), center = TRUE, scale. = FALSE)

percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  analysis_group = colData(vsd)$analysis_group
)

pca_colors <- ifelse(pca_df$analysis_group == "Tumor", "firebrick", "steelblue")

png(filename = pca_plot_path, width = 2400, height = 1800, res = 300)
plot(
  pca_df$PC1,
  pca_df$PC2,
  col = pca_colors,
  pch = 16,
  cex = 0.8,
  xlab = paste0("PC1: ", percent_var[1], "% variance"),
  ylab = paste0("PC2: ", percent_var[2], "% variance"),
  main = "TCGA-COAD PCA: Tumor vs Normal"
)
legend(
  "topright",
  legend = c("Normal", "Tumor"),
  col = c("steelblue", "firebrick"),
  pch = 16,
  title = "Group"
)
dev.off()

# ---------------------------
# 6. Volcano plot
# ---------------------------

volcano_df <- res_df %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "Upregulated in tumor",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated in tumor",
      TRUE ~ "Not significant"
    ),
    neg_log10_padj = -log10(padj),
    abs_log2fc = abs(log2FoldChange)
  )

# Replace infinite values if any
volcano_df$neg_log10_padj[is.infinite(volcano_df$neg_log10_padj)] <- NA_real_

# Keep only rows that can be plotted
volcano_plot_df <- volcano_df %>%
  filter(!is.na(log2FoldChange), !is.na(neg_log10_padj))

# Select more genes to label
top_up <- volcano_plot_df %>%
  filter(significance == "Upregulated in tumor") %>%
  arrange(padj, desc(abs_log2fc)) %>%
  slice_head(n = 6)

top_down <- volcano_plot_df %>%
  filter(significance == "Downregulated in tumor") %>%
  arrange(padj, desc(abs_log2fc)) %>%
  slice_head(n = 6)

# Helper function to reduce overlap by spacing labels vertically
adjust_label_positions <- function(df, min_gap = 8) {
  if (nrow(df) <= 1) {
    df$label_y <- df$neg_log10_padj
    return(df)
  }

  df <- df[order(df$neg_log10_padj), ]
  df$label_y <- df$neg_log10_padj

  for (i in 2:nrow(df)) {
    if ((df$label_y[i] - df$label_y[i - 1]) < min_gap) {
      df$label_y[i] <- df$label_y[i - 1] + min_gap
    }
  }

  df
}

top_up <- adjust_label_positions(top_up, min_gap = 8)
top_down <- adjust_label_positions(top_down, min_gap = 8)

volcano_labels <- bind_rows(top_up, top_down)

# Plot non-significant points first, then significant ones
volcano_nonsig <- volcano_plot_df %>% filter(significance == "Not significant")
volcano_down <- volcano_plot_df %>% filter(significance == "Downregulated in tumor")
volcano_up <- volcano_plot_df %>% filter(significance == "Upregulated in tumor")

# Define plotting ranges with a small buffer
x_range <- range(volcano_plot_df$log2FoldChange, na.rm = TRUE)
y_range <- range(volcano_plot_df$neg_log10_padj, na.rm = TRUE)
label_y_max <- max(volcano_labels$label_y, na.rm = TRUE)

x_buffer <- 0.5
y_buffer <- 5

png(filename = volcano_plot_path, width = 2400, height = 1800, res = 300)

plot(
  volcano_nonsig$log2FoldChange,
  volcano_nonsig$neg_log10_padj,
  col = "grey75",
  pch = 16,
  cex = 0.5,
  xlab = "Log2 fold change (Tumor vs Normal)",
  ylab = expression(-log[10]("adjusted p-value")),
  main = "TCGA-COAD Volcano Plot",
  xlim = c(x_range[1] - x_buffer, x_range[2] + x_buffer),
  ylim = c(0, max(y_range[2], label_y_max) + y_buffer)
)

points(
  volcano_down$log2FoldChange,
  volcano_down$neg_log10_padj,
  col = "steelblue",
  pch = 16,
  cex = 0.6
)

points(
  volcano_up$log2FoldChange,
  volcano_up$neg_log10_padj,
  col = "firebrick",
  pch = 16,
  cex = 0.6
)

abline(v = c(-1, 1), lty = 2)
abline(h = -log10(0.05), lty = 2)

# Draw connector lines from points to labels
segments(
  x0 = volcano_labels$log2FoldChange,
  y0 = volcano_labels$neg_log10_padj,
  x1 = volcano_labels$log2FoldChange,
  y1 = volcano_labels$label_y,
  col = "black",
  lty = 1
)

# Add labels
with(
  volcano_labels,
  text(
    x = log2FoldChange,
    y = label_y,
    labels = plot_label,
    pos = ifelse(log2FoldChange > 0, 4, 2),
    cex = 0.75,
    offset = 0.35,
    xpd = NA
  )
)

legend(
  "topleft",
  inset = 0.02,
  legend = c("Upregulated in tumor", "Downregulated in tumor", "Not significant"),
  col = c("firebrick", "steelblue", "grey75"),
  pch = 16,
  bty = "n"
)

dev.off()

# ---------------------------
# 7. MA plot
# ---------------------------

png(filename = ma_plot_path, width = 2400, height = 1800, res = 300)
plotMA(
  res,
  ylim = c(-5, 5),
  main = "TCGA-COAD MA Plot"
)
dev.off()

# ---------------------------
# 8. Heatmap of top 30 genes
# ---------------------------

top_genes_df <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 30)

write.csv(top_genes_df, file = top_genes_csv_path, row.names = FALSE)

top_gene_ids <- top_genes_df$gene_id
top_gene_labels <- top_genes_df$plot_label

heatmap_matrix <- assay(vsd)[top_gene_ids, , drop = FALSE]

# Row-scale expression values for visualization
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
rownames(heatmap_matrix_scaled) <- make.unique(top_gene_labels)

annotation_col <- as.data.frame(colData(vsd)[, "analysis_group", drop = FALSE])
annotation_col$analysis_group <- as.character(annotation_col$analysis_group)
rownames(annotation_col) <- colnames(vsd)

pheatmap(
  mat = heatmap_matrix_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Top 30 Differentially Expressed Genes",
  filename = heatmap_plot_path,
  width = 10,
  height = 8
)

# ---------------------------
# 9. Final summary
# ---------------------------

cat("Visualization script completed successfully.\n")
cat("Saved files:\n")
cat("- ", pca_plot_path, "\n", sep = "")
cat("- ", volcano_plot_path, "\n", sep = "")
cat("- ", ma_plot_path, "\n", sep = "")
cat("- ", heatmap_plot_path, "\n", sep = "")
cat("- ", top_genes_csv_path, "\n", sep = "")