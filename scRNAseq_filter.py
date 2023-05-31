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
    parser.add_argument('p_value', help='p-value threshold to define outliers for QC filtering removal', type=float)
   
    # Parse args
    args = parser.parse_args()

    # directory that contains features, barcode, and matrix files
    data_dir = args.dir 
    p_value = args.p_value
    
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

    # get cutoff value for p < p_value
    print()
    print(f'n_gene_by_counts cutoff value (p < {p_value})', n_genes_by_counts.find_cutoff(n_genes_by_counts_dist, p_value))
    print(f'total_counts cutoff value (p < {p_value})', total_counts.find_cutoff(total_counts_dist, p_value))
    print(f'pct_counts_mt cutoff value (p < {p_value})', pct_counts_mt.find_cutoff(pct_counts_mt_dist, p_value))
    return

if __name__ == "__main__":
    main()
    
