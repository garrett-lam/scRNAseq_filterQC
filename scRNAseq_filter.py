import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import preprocess_data
import n_genes_by_counts
import total_counts
import pct_counts_mt

# MAIN DRIVER
def main():
    # Create arg parser obj
    parser = argparse.ArgumentParser(
        prog="scRNAseq_filter",
        description="Command-line script to perform filtering of scRNA-seq data"
    )
    # Positional Args
    parser.add_argument('dir', help='Directory containing scRNA-seq barcodes, features, matrix files', type=str)
    # Parse args
    args = parser.parse_args()

    data_dir = args.dir 
    adata_obj = preprocess_data.preprocess_data(data_dir)
    n_genes_by_counts.cutoff_005(adata_obj)
    total_counts.cutoff_005(adata_obj)
    pct_counts_mt.cutoff_005(adata_obj)

    return

if __name__ == "__main__":
    main()
    
