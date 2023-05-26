# scRNAseq_filterQC

**Motivation/Description of Project:** As demonstrated in Lab 6, single-cell RNA-seq (scRNA-seq) data quality control has mostly been done qualitatively. Thresholds to filter out outliers are highly subjective and selected at the user's discretion via data visualizations. As a result, we wanted to take a more quantitative approach to standardize this process of outlier-removal. This project is focused on developing a python tool to streamline the process of quality control filtering of scRNA-seq data. Specifically, we will apply statistical models and/or transformations to identify discrete thresholds for outlier removal (high mitochondrial read content, genes by counts, total counts). 

## Install Instructions & Environment Setup
There are 2 ways to set up the environment: 
- Create environment using `conda env create -f cse185.yaml` then `conda activate cse185` to activate the environment 
- Create environment using `micromamba env create --name cse185 --file cse185_micromamba.yaml` then `micromamba activate cse185` to activate the environment.

## Manual Preprocessing Steps:
**IMPORTANT:** Ensure barcodes, features, matrix files are formatted in the following way:

`[filename]_barcodes.tsv`
- **Column 1:** Barcodes

`[filename]_features.tsv`
- **Column 1:** Gene IDs
- **Column 2:** Gene Symbols
- **Column 3:** The string: 'Gene Expression'

`[filename]_matrix.mtx`
First line contains dimensions of the matrix and total number of non-zero entries in the format: rows, columns, non-zero entries
Following lines after the first two lines must be numeric values
- **Column 1:** Row index of matrix (gene/feature)
- **Column 2:** Column index of matrix (cell)
- **Column 3:** Value of matrix entry (expression level of gene in cell)

## Basic Usage:
Ensure barcodes, features, matrix files are put into `data_dir`:
`python scRNAseq_filter.py [data_dir]`

## How to run test example (`pbmc_test`):
`python scRNAseq_filter.py pbmc_test/`

## Plan
 1. Fit count QC metrics (`n_genes_by_count`, `total_counts`) to the Poisson Distribution
 2. Fit continuous QC metrics (`pct_counts_mt`) to the Weibull Distribution
 3. For the 3 above QC metrics, do the following to handle outliers:
    - For data points with a p-value < 0.05 or p-value < 0.01, winsorize them or remove them
 4. Perform clustering of the cells using the Leiden algorithm
 5. Assess clustering results using t-SNE and UMAP visualizations and objective measures

