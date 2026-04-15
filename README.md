# TCGA-COAD Bulk RNA-seq Differential Expression Analysis

A reproducible downstream bulk RNA-seq analysis of TCGA colon adenocarcinoma (TCGA-COAD) to identify differentially expressed genes between tumor and normal samples using DESeq2 in R.

## Overview
This project explores transcriptomic differences between tumor and normal samples in TCGA-COAD using a reproducible differential expression workflow. The analysis is designed to demonstrate practical skills in human cancer transcriptomics, RNA-seq data handling, differential expression analysis, visualization, and biological interpretation.

## Project objective
The main objective of this project is to identify genes that are differentially expressed between colon adenocarcinoma tumor samples and normal colon samples, and to summarize the resulting transcriptomic patterns through clear and reproducible analyses.

## Dataset
- **Source:** TCGA-COAD
- **Data type:** Bulk RNA-seq gene expression data
- **Comparison:** Tumor vs normal samples

## Tools and packages
- R
- TCGAbiolinks
- DESeq2
- ggplot2
- pheatmap

## Planned analysis workflow
1. Retrieve TCGA-COAD RNA-seq count data and associated sample metadata
2. Prepare and clean metadata for downstream analysis
3. Perform differential expression analysis with DESeq2
4. Generate exploratory and result-focused visualizations
5. Interpret key transcriptomic findings in a biological context

## Expected outputs
This project will include:
- metadata summary tables
- quality control and exploratory plots
- differential expression results
- visualization of significant transcriptomic patterns
- a reproducible and clearly documented analysis structure

## Repository structure
- `data/` — input and processed data files
- `scripts/` — analysis scripts
- `results/` — statistical outputs and result tables
- `figures/` — plots and visual summaries

## Status
Project setup in progress.
