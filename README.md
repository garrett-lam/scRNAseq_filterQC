# scRNAseq_filterQC

## Environment Setup
There are 2 ways to set up the environment: `conda env create -f cse185.yaml` or `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Depending on the method used, use `conda activate cse185` or `micromamba activate cse185` to activate the environment.

## Notes
 - Call `gunzip anndata_with_QCMetrics.h5ad.gz` to use the anndata object.

## Plan
 1. Fit count QC metrics (`n_genes_by_count`, `total_counts`) to the Poisson Distribution
 2. Fit continuous QC metrics (``pct_counts_mt`) to the Weibull Distribution
 3. For the 3 above QC metrics, do the following to handle outliers:
    - For data points with a p-value < 0.05 or p-value < 0.01, winsorize them or remove them
 4. Perform clustering of the cells using the Leiden algorithm
 5. Assess clustering results using t-SNE and UMAP visualizations and objective measures
