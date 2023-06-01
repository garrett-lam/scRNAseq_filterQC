import argparse
import preprocess_data
import n_genes_by_counts
import total_counts
import pct_counts_mt
import scanpy as sc
import scanpy.external as sce
import warnings

# MAIN DRIVER
def main():
    warnings.filterwarnings("ignore")
    warnings.filterwarnings("ignore", category=DeprecationWarning, module=".*mkl.*")

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

    n_genes_by_counts_cutoff = n_genes_by_counts.find_cutoff(n_genes_by_counts_dist, p_value)
    total_counts_cutoff = total_counts.find_cutoff(total_counts_dist, p_value)
    pct_counts_mt_cutoff = pct_counts_mt.find_cutoff(pct_counts_mt_dist, p_value)

    print(f'n_gene_by_counts cutoff value (p < {p_value})', n_genes_by_counts_cutoff)
    print(f'total_counts cutoff value (p < {p_value})', total_counts_cutoff)
    print(f'pct_counts_mt cutoff value (p < {p_value})', pct_counts_mt_cutoff)

    #filter the QC metrics 
    adata_obj = adata_obj[adata_obj.obs.n_genes_by_counts < n_genes_by_counts_cutoff, :] # SET TO NEW VARIABLE?
    adata_obj = adata_obj[adata_obj.obs.total_counts < total_counts_cutoff, :]
    adata_obj = adata_obj[adata_obj.obs.pct_counts_mt < pct_counts_mt_cutoff, :]

    print(f'Number of Genes (after filtering): {adata_obj.n_vars}')
    print(f'Number of Cells (after filtering): {adata_obj.n_obs}')

    sc.pp.normalize_per_cell(adata_obj, counts_per_cell_after=1e4) # normalize to 10,000 reads/cell
    sc.pp.log1p(adata_obj) # log transform

    # get top 500 variable genes
    sc.pp.highly_variable_genes(adata_obj, batch_key="dataset", n_top_genes=500)

    print(adata_obj.var[adata_obj.var['highly_variable']].sort_values(by='dispersions_norm', ascending=False)) # OPTIONAL

    ''' (FROM LAB6 )
    # For the analyses below, we recommend making a new anndata object, which contains only:
    - Highly variable genes
    - Genes in the set of cell-type specific marker genes used in the paper (see below). We will manually add these back, since we want to analyze them even if they didn't make the cut for being most differentially expressed.
    '''
    # Run PCA
    sc.pp.pca(adata_obj, n_comps=20)
    #sc.pl.pca(adata_obj, color="dataset", title='PCA plot before removing batch effects')

    # Remove batch effects
    sce.pp.harmony_integrate(adata_obj, 'dataset', theta=2, nclust=50,  max_iter_harmony = 10,  max_iter_kmeans=10)

    # Reset the original PCs to those computed by Harmony
    adata_obj.obsm['X_pca'] = adata_obj.obsm['X_pca_harmony']

    # Perform clustering
    sc.pp.neighbors(adata_obj) # computes neighborhood graphs. Needed to run clustering.
    sc.tl.leiden(adata_obj) # clusters cells based on expression profiles. This is needed to color cells by cluster. # ERROR HERE

    # UMAP 
    sc.tl.umap(adata_obj) # compute UMAP embedding
    sc.pl.umap(adata_obj, color="leiden") # plot UMAP, coloring cells by cluster
    sc.pl.umap(adata_obj, color="dataset")  # plot UMAP, coloring cells by dataset

    # we can include tSNE if you want but I think UMAP is better lol

    return

if __name__ == "__main__":
    main()
