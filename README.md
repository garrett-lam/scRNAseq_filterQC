# scRNAseqFilterQC

**Motivation/Description of Project:** As demonstrated in Lab 6, single-cell RNA-seq (scRNA-seq) data quality control has mostly been done qualitatively. Thresholds to filter out outliers are highly subjective and selected at the user's discretion via data visualizations. As a result, we wanted to take a more quantitative approach to standardize this process of outlier-removal. This project is focused on developing a python tool to streamline the process of quality control filtering of scRNA-seq data. Specifically, we will apply statistical models to identify discrete thresholds for outlier removal. We will start by focusing on the following QC metrics: percent mitochondrial read content, genes by counts, and total counts. 

## Install Instructions & Environment Setup
First, clone the repo by calling `git clone https://github.com/garrett-lam/scRNAseq_filterQC.git`. Then, change into the directory using `cd scRNAseq_filterQC`.

Then, there are 2 ways to set up the environment (you only need to do one of these): 
 - __If you use/prefer conda__: Create environment using `conda env create -f cse185.yaml`. Then call `conda activate cse185` to activate the environment.
 - __If you use/prefer micromamba__: Create environment using `micromamba env create --name cse185 --file cse185_micromamba.yaml`. Then call `micromamba activate cse185` to activate the environment.

___Note for Datahub Users___:
 - If you get `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.` when calling `conda activate cse185`, follow the onscreen directions by running `conda init bash` (that is the shell used on Datahub) and closing/restarting your shell. Then, instead of calling `conda activate cse185`, call `source activate cse185` to activate the environment.

___General Notes___:
 - The environment creation step may take a while. This is normal.
 - If you are prompted with a message to update conda `==> WARNING: A newer version of conda exists. <==`, you may safely ignore it.
 - If you are prompted with a warning, `Warning: you have pip-installed dependences in your environment file, but you do not list pip itself as one of your conda dependencies...`, you may safely ignore it.

## Manual Preprocessing Steps:
**IMPORTANT:** Ensure barcodes, features, matrix files are formatted in the following way:

`[filename]_barcodes.tsv`
 - **Column 1:** Barcodes

`[filename]_features.tsv`
 - **Column 1:** Gene IDs
 - **Column 2:** Gene Symbols
 - **Column 3:** The string: 'Gene Expression'

`[filename]_matrix.mtx`
- First line contains `%%MatrixMarket matrix coordinate integer general`
- Second line contains dimensions of the matrix and total number of non-zero entries in the format: rows, columns, non-zero entries

Following lines after the first two lines must be numeric values
 - **Column 1:** Row index of matrix (gene/feature)
 - **Column 2:** Column index of matrix (cell)
 - **Column 3:** Value of matrix entry (expression level of gene in cell)

## Basic Usage:
Call `python scRNAseqFilterQC.py data_dir/ [other options]`
- `data_dir/` is the directory that contains barcodes, features, and matrix files. 

## scRNAseqFilterQC options
- `-n'` or `--n_genes_by_counts_p_value`: p-value threshold to define outliers for QC filtering removal by `n_genes_by_counts`
- `-t` or `--total_counts_p_value`: p-value threshold to define outliers for QC filtering removal by `total_counts`
- `-p` or `--pct_counts_mt_p_value`: p-value threshold to define outliers for QC filtering removal by `pct_counts_mt`
- `-g` or `--marker_genes`: Cell-type specific marker genes of interest

## How to run the workflow using the dataset:
**Lab 6 Data**
```
python scRNAseqFilterQC.py counts/ -n 0.01 -t 0.01 -p 0.125 \
    -g GCG TTR IAPP GHRL PPY COL3A1 CPA1 CLPS REG1A CTRB1 CTRB2 PRSS2 CPA2 KRT19 INS SST CELA3A VTCN1
```

**Real World Data (pbmc3k)**
```
python scRNAseqFilterQC.py pbmc_test/ -n 0.01 -t 0.01 -p 0.125 \
    -g IL7R CD79A MS4A1 CD8A CD8B LYZ CD14 LGALS3 S100A8 GNLY NKG7 KLRB1 FCGR3A MS4A7 FCER1A CST3 PPBP
```

