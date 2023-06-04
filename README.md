# scRNAseqFilterQC

__Motivation/Description of Project:__ As demonstrated in Lab 6, single-cell RNA-seq (scRNA-seq) data quality control has mostly been done qualitatively. Thresholds to filter out outliers are highly subjective and selected at the user's discretion via data visualizations. As a result, we wanted to take a more quantitative approach to standardize this process of outlier-removal. This project is focused on developing a python tool to streamline the process of quality control filtering of scRNA-seq data. Specifically, we will apply statistical models to identify discrete thresholds for outlier removal. We will start by focusing on the following QC metrics: percent mitochondrial read content, genes by counts, and total counts. 

## Install Instructions & Environment Setup
First, clone the repo by calling `git clone https://github.com/garrett-lam/scRNAseq_filterQC.git`. Then, change into the directory using `cd scRNAseq_filterQC`.

Then, there are 2 ways to set up the environment (you only need to do one of these): 
 - __If you use/prefer conda__: Create environment using `conda env create -f cse185.yaml`. Then call `conda activate cse185` to activate the environment.
 - __If you use/prefer micromamba__: Create environment using `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Then call `micromamba activate cse185` to activate the environment.

___Note for Datahub Users___:
 - If you get `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.` when calling `conda activate cse185`, follow the onscreen directions by running `conda init bash` (that is the shell used on Datahub). After this, run the `exit` command, and then open a new Terminal. In this new Terminal, call  `source activate cse185` instead of `conda activate cse185` to activate the environment. ___If this happens to you, always use `source activate cse185` instead of `conda activate cse185`___.

___General Notes___:
 - The environment creation step may take a while. This is normal.
 - If you are prompted with a message to update conda `==> WARNING: A newer version of conda exists. <==`, you may safely ignore it.
 - If you are prompted with a warning, `Warning: you have pip-installed dependences in your environment file, but you do not list pip itself as one of your conda dependencies...`, you may safely ignore it.
 - If you have weird errors with the `conda` environment (e.g., running into segfaults, code stalling, etc.), set up the `micromamba` environment and try running the code from there instead.

## Manual Preprocessing Steps:
__IMPORTANT:__ Ensure barcodes, features, matrix files are formatted in the following way:

`[filename]_barcodes.tsv`
 - __Column 1:__ Barcodes

`[filename]_features.tsv`
 - __Column 1:__ Gene IDs
 - __Column 2:__ Gene Symbols
 - __Column 3:__ The string: 'Gene Expression'

`[filename]_matrix.mtx`
- First line contains `%%MatrixMarket matrix coordinate integer general`
- Second line contains dimensions of the matrix and total number of non-zero entries in the format: rows, columns, non-zero entries

Following lines after the first two lines must be numeric values
 - __Column 1:__ Row index of matrix (gene/feature)
 - __Column 2:__ Column index of matrix (cell)
 - __Column 3:__ Value of matrix entry (expression level of gene in cell)

## Basic Usage:
Call `python scRNAseqFilterQC.py data_dir/ [other options]`
- `data_dir/` is the directory that contains barcodes, features, and matrix files. Note that when passing `data_dir`, you must include the `/` after the directory name.

## scRNAseqFilterQC options
- `-n` or `--n_genes_by_counts_p_value`: p-value threshold to define outliers for QC filtering removal by `n_genes_by_counts`
- `-t` or `--total_counts_p_value`: p-value threshold to define outliers for QC filtering removal by `total_counts`
- `-p` or `--pct_counts_mt_p_value`: p-value threshold to define outliers for QC filtering removal by `pct_counts_mt`
- `-g` or `--marker_genes`: Cell-type specific marker genes of interest

## Output Details
The script will create an output folder named as follows:
```
<data_dir>_<n_genes_by_counts_p_value>_<total_counts_p_value>_<pct_counts_md_p_value>
```

_Note that `<data_dir>` in the output folder name does not include the `/`._

Within this folder, you can fine the following files:
 - `figures/best_dist_n_genes_by_counts.png`: top 5 distributions fitted to `n_genes_by_counts`, sorted low to high by MSE in the legend. The first distribution listed represents the final distribution the data were fit to.
 - `figures/best_dist_pct_counts_mt.png`: top 5 distributions fitted to `pct_counts_mt`, sorted low to high by MSE in the legend. The first distribution listed represents the final distribution the data were fit to.
 - `figures/best_dist_total_counts.png`: top 5 distributions fitted to `total_counts`, sorted low to high by MSE in the legend. The first distribution listed represents the final distribution the data were fit to.
 - `figures/tsne_clusters.png`: t-SNE plot colored by assigned clusters
 - `figures/tsne_datasets.png`: t-SNE plot colored by dataset
 - `figures/tsne_marker_genes.png`: t-SNE plot colored by selected marker genes
 - `figures/umap_clusters.png`: UMAP plot colored by assigned clusters
 - `figures/umap_datasets.png`: UMAP plot colored by dataset
 - `figures/umap_marker_genes.png`: UMAP plot colored by selected marker genes

## How to run the workflow using the dataset:
__Lab 6 Data__
```
python scRNAseqFilterQC.py counts/ -n 0.01 -t 0.01 -p 0.125 \
    -g GCG TTR IAPP GHRL PPY COL3A1 CPA1 CLPS REG1A CTRB1 CTRB2 PRSS2 CPA2 KRT19 INS SST CELA3A VTCN1
```

__Real World Data (pbmc3k)__
```
python scRNAseqFilterQC.py pbmc_test/ -n 0.01 -t 0.01 -p 0.01 \
    -g IL7R CD79A MS4A1 CD8A CD8B LYZ CD14 LGALS3 S100A8 GNLY NKG7 KLRB1 FCGR3A MS4A7 FCER1A CST3 PPBP
```

