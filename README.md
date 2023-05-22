# scRNAseq_filterQC

## Environment Setup
There are 2 ways to set up the environment: `conda env create -f cse185.yaml` or `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Depending on the method used, use `conda activate cse185` or `micromamba activate cse185` to activate the environment.

## Manual Preprocessing Steps:
**IMPORTANT:** Ensure barcodes, features, matrix files are formatted in the following way:

`[filename]_barcodes.tsv`
- Column 1: Barcodes

`[filename]_features.tsv`
- Column 1: Gene IDs
- Column 2: Gene Symbols
- Column 3: 'Gene Expression'

`[filename]_matrix.mtx`
- First line in file must be: %%MatrixMarket matrix coordinate real general 
- Second line contains dimensions of the matrix and total number of non-zero entries in the format: rows columns non-zero entries

Following lines after the first two lines must be numeric values
- Column 1: Row index of matrix (gene/feature)
- Column 2: Column index of matrix (cell)
- Column 3: Value of matrix entry (expression level of gene in cell)

## How to run test example (`pbmc_test`):
`python scRNAseq_filter.py pbmc_test/`

## Plan
 1. Fit count QC metrics (`n_genes_by_count`, `total_counts`) to the Poisson Distribution
 2. Fit continuous QC metrics (``pct_counts_mt`) to the Weibull Distribution
 3. For the 3 above QC metrics, do the following to handle outliers:
    - For data points with a p-value < 0.05 or p-value < 0.01, winsorize them or remove them
 4. Perform clustering of the cells using the Leiden algorithm
 5. Assess clustering results using t-SNE and UMAP visualizations and objective measures



