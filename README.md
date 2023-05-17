# scRNAseq_filterQC

## Environment Setup
There are 2 ways to set up the environment: `conda env create -f cse185.yaml` or `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Depending on the method used, use `conda activate cse185` or `micromamba activate cse185` to activate the environment.

## Notes
 - Call `gunzip anndata_with_QCMetrics.h5ad.gz` to use the anndata object.
