# scRNAseq_filterQC

## Environment Setup
There are 2 ways to set up the environment: `conda env create -f cse185.yaml` or `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Depending on the method used, use `conda activate cse185` or `micromamba activate cse185` to activate the environment.

## Notes
 - Call `gunzip anndata_with_QCMetrics.h5ad.gz` to use the anndata object.


## Manual Preprocessing Steps:
- Ensure barcodes, features, matrix files are formatted in the following way:

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