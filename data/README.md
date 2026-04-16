# Data folder

This folder stores downloaded and generated data files used in the TCGA-COAD bulk RNA-seq differential expression workflow.

Large data files are not tracked in GitHub in order to keep the repository lightweight and reproducible.

## How to generate the data

Run the following script from the project root:

```bash
Rscript scripts/01_download_tcga_coad_data.R