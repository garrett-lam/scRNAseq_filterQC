import argparse
import preprocess_data
import n_genes_by_counts
import total_counts
import pct_counts_mt
import scanpy as sc
import scanpy.external as sce
import warnings
import subprocess
import os

# MAIN DRIVER
def main():
    
    warnings.filterwarnings("ignore") # Filter out warnings

    # Create arg parser obj
    parser = argparse.ArgumentParser(
        prog="scRNAseq_filter",
        description="Command-line script to perform filtering of scRNA-seq data"
    )
    
    # Positional Args
    parser.add_argument('dir', help='Directory containing scRNA-seq barcodes, features, matrix files', type=str)
    parser.add_argument('-n', '--n_genes_by_counts_p_value', help='n_genes_by_counts p-value threshold to define outliers for QC filtering removal', type=float)
    parser.add_argument('-t', '--total_counts_p_value', help='total_counts p-value threshold to define outliers for QC filtering removal', type=float)
    parser.add_argument('-p', '--pct_counts_mt_p_value', help='pct_counts_mt p-value threshold to define outliers for QC filtering removal', type=float)
    parser.add_argument('-g', '--marker_genes', help='Cell-type specific marker genes of interest', nargs='*', type=str, default=[]) # 0 or more occurrences

    # Parse args
    args = parser.parse_args()

    data_dir = args.dir # Directory that contains features, barcode, and matrix files
    n_genes_by_counts_p_value = args.n_genes_by_counts_p_value
    total_counts_p_value = args.total_counts_p_value
    pct_counts_mt_p_value = args.pct_counts_mt_p_value
    marker_genes = args.marker_genes # list of strings
    
    # Preprocess data
    adata_obj = preprocess_data.preprocess_data(data_dir)

    out_dir = data_dir[:-1] + '_' + str(n_genes_by_counts_p_value) + '_' + str(total_counts_p_value) + '_' + str(pct_counts_mt_p_value)
    subprocess.run(["mkdir", "-p", out_dir])
    os.chdir(out_dir)

    # Find the best distribution for the three catagories
    subprocess.run(["mkdir", "-p", "figures"])
    n_genes_by_counts_dist = n_genes_by_counts.find_distributions(adata_obj)
    total_counts_dist = total_counts.find_distributions(adata_obj)
    pct_counts_mt_dist = pct_counts_mt.find_distributions(adata_obj)

    # Return the best distribution for the three catagories
    print(f'n_genes_by_counts: {n_genes_by_counts_dist}')
    print(f'total_counts: {total_counts_dist}')
    print(f'pct_counts_mt: {pct_counts_mt_dist}')

    # Get cutoff value for p < p_value
    print()

    n_genes_by_counts_cutoff = n_genes_by_counts.find_cutoff(n_genes_by_counts_dist, n_genes_by_counts_p_value)
    total_counts_cutoff = total_counts.find_cutoff(total_counts_dist, total_counts_p_value)
    pct_counts_mt_cutoff = pct_counts_mt.find_cutoff(pct_counts_mt_dist, pct_counts_mt_p_value)

    print(f'n_gene_by_counts cutoff value (p < {n_genes_by_counts_p_value})', n_genes_by_counts_cutoff)
    print(f'total_counts cutoff value (p < {total_counts_p_value})', total_counts_cutoff)
    print(f'pct_counts_mt cutoff value (p < {pct_counts_mt_p_value})', pct_counts_mt_cutoff)

    #filter the QC metrics 
    adata_obj = adata_obj[adata_obj.obs.n_genes_by_counts < n_genes_by_counts_cutoff, :] 
    adata_obj = adata_obj[adata_obj.obs.total_counts < total_counts_cutoff, :]
    adata_obj = adata_obj[adata_obj.obs.pct_counts_mt < pct_counts_mt_cutoff, :]

    print(f'Number of Genes (after filtering): {adata_obj.n_vars}')
    print(f'Number of Cells (after filtering): {adata_obj.n_obs}')

    sc.pp.normalize_per_cell(adata_obj, counts_per_cell_after=1e4) # normalize to 10,000 reads/cell
    sc.pp.log1p(adata_obj) # log transform

    # get top 500 variable genes
    sc.pp.highly_variable_genes(adata_obj, batch_key="dataset", n_top_genes=500)

    print(adata_obj.var[adata_obj.var['highly_variable']].sort_values(by='dispersions_norm', ascending=False)) 
    
    if marker_genes:
        adata_filt = adata_obj[:, (adata_obj.var.index.isin(marker_genes) | adata_obj.var["highly_variable"])]
    else:
        adata_filt = adata_obj[:, (adata_obj.var["highly_variable"])]

    # Run PCA
    sc.pp.pca(adata_filt, n_comps=20)

    # Remove batch effects
    sce.pp.harmony_integrate(adata_filt, 'dataset', theta=2, nclust=50,  max_iter_harmony = 10,  max_iter_kmeans=10)

    # Reset the original PCs to those computed by Harmony
    adata_filt.obsm['X_pca'] = adata_filt.obsm['X_pca_harmony']

    # Perform clustering
    sc.pp.neighbors(adata_filt) # computes neighborhood graphs. Needed to run clustering.
    sc.tl.leiden(adata_filt) # clusters cells based on expression profiles. This is needed to color cells by cluster. # ERROR HERE

    # UMAP 
    sc.tl.umap(adata_filt) # compute UMAP embedding
    sc.pl.umap(adata_filt, color="leiden", save='_clusters.png') # plot UMAP, coloring cells by cluster
    sc.pl.umap(adata_filt, color="dataset", save='_datasets.png')  # plot UMAP, coloring cells by dataset
    sc.pl.umap(adata_filt, color=marker_genes, color_map="Reds", save='_marker_genes.png')

    # tSNE
    sc.tl.tsne(adata_filt)
    sc.pl.tsne(adata_filt, color=['leiden'], legend_loc='on data', legend_fontsize=10, alpha=0.8, size=20, save='_clusters.png')
    sc.pl.tsne(adata_filt, color=['dataset'], legend_loc='on data', legend_fontsize=10, alpha=0.8, size=20, save='_datasets.png')
    sc.pl.tsne(adata_filt, color=marker_genes, color_map="Reds", save='_marker_genes.png')
    
    return

if __name__ == "__main__":
    main()
