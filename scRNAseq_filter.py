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

    # directory that contains features, barcode, and matrix files
    data_dir = args.dir 
    
    # preprocess data
    adata_obj = preprocess_data.preprocess_data(data_dir)

    # find the best distribution for the three catagories
    n_genes_by_counts_dist = n_genes_by_counts.find_distributions(adata_obj)
    total_counts_dist = total_counts.find_distributions(adata_obj)
    pct_counts_mt_dist = pct_counts_mt.find_distributions(adata_obj)

    # return the best distribution for the three catagories
    print(f'n_genes_by_counts: {n_genes_by_counts_dist}')
    print(f'total_counts: {total_counts_dist}')
    print(f'pct_counts_mt: {pct_counts_mt_dist}')
    return

if __name__ == "__main__":
    main()
    
