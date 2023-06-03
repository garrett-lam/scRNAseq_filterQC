
import sys
import os
import scanpy as sc
import harmonypy
import leidenalg
import anndata as ad
# Import the "external" library
import scanpy.external as sce



sys.path.append(os.environ["HOME"]+"/.local/lib/python3.9/site-packages")
sc.logging.print_versions()
DATADIR=os.environ["HOME"]+"/public/lab6"
dataset = sc.read_10x_mtx(DATADIR, prefix="GSM5114461_S6_A11_", cache=True)
DATADIR=os.environ["HOME"]+"/public/lab6"
dsets = ["GSM5114461_S6_A11", "GSM5114464_S7_D20", "GSM5114474_M3_E7"]
adatas = {}
for ds in dsets:
    print(ds)
    adatas[ds] = sc.read_10x_mtx(DATADIR, prefix=ds+"_", cache=True)
combined = ad.concat(adatas, label="dataset")
combined.obs_names_make_unique()
total_cells = 0
for ds, adata in adatas.items():
    n_genes = adata.shape[1]
    n_cells = adata.shape[0]
    print(f"Dataset {ds} has {n_cells} cells and {n_genes} genes.")
    total_cells += n_cells

print(f"The combined dataset has {total_cells} cells and {n_genes} genes.")
sc.pp.filter_cells(combined, min_genes=200)
sc.pp.filter_cells(combined, min_counts=1000)

# Filter genes
sc.pp.filter_genes(combined, min_cells=5)
sc.pp.filter_genes(combined, min_counts=15)

# Print the number of cells and genes after filtering
print(f"Number of cells after filtering: {combined.shape[0]}")
print(f"Number of genes after filtering: {combined.shape[1]}")
sc.pl.highest_expr_genes(combined, n_top=20)

combined.var['mt'] = combined.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(combined, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(combined, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(combined, x='total_counts', y='n_genes_by_counts')
combined_filt = combined[combined.obs.pct_counts_mt < 25, :]
combined_filt = combined_filt[combined_filt.obs.total_counts < 70000, :]
print(f"Number of cells after filtering: {combined_filt.shape[0]}")
print(f"Number of genes after filtering: {combined_filt.shape[1]}")

sc.pp.normalize_per_cell(combined_filt, counts_per_cell_after=1e4) # normalize to 10,000 reads/cell
sc.pp.log1p(combined_filt)

# Identify highly variable genes
sc.pp.highly_variable_genes(combined_filt, batch_key="dataset", n_top_genes=500)
combined_filt.var[combined_filt.var['highly_variable']].sort_values(by='dispersions_norm', ascending=False)

genes = ["GCG", "TTR",  "IAPP",  "GHRL", "PPY", "COL3A1",
    "CPA1", "CLPS", "REG1A", "CTRB1", "CTRB2", "PRSS2", "CPA2", "KRT19", "INS","SST","CELA3A", "VTCN1"]

combined_var = combined_filt[:, (combined_filt.var.index.isin(genes) | combined_filt.var["highly_variable"])]

sc.pp.pca(combined_var, n_comps = 20)
sc.pl.pca(combined_var, color="dataset")



# Run harmony using suggested params from the paper
sce.pp.harmony_integrate(combined_var, 'dataset', theta=2, nclust=50,  max_iter_harmony = 10,  max_iter_kmeans=10)

# Reset the original PCs to those computed by Harmony
combined_var.obsm['X_pca'] = combined_var.obsm['X_pca_harmony']
sc.pl.pca(combined_var, color="dataset")
sc.pp.neighbors(combined_var) # computes neighborhood graphs. Needed to run clustering.
sc.tl.leiden(combined_var) # clusters cells based on expression profiles. This is needed to color cells by cluster.
sc.tl.umap(combined_var) # compute UMAP embedding
sc.pl.umap(combined_var, color="leiden") # make the UMAP plot, coloring cells by cluster
sc.pl.umap(combined_var, color="dataset") # make the UMAP plot, coloring cells by cluster
sc.tl.tsne(combined_var)
sc.pl.tsne(combined_var, color=['leiden'], legend_loc='on data', legend_fontsize=10, alpha=0.8, size=20)
sc.pl.tsne(combined_var, color=['dataset'], legend_loc='on data', legend_fontsize=10, alpha=0.8, size=20)
sc.pl.umap(combined_var, color=genes, color_map="Reds")
sc.pl.tsne(combined_var, color=genes, color_map="Reds")
#sc.pl.heatmap(combined_var, genes, groupby='leiden', dendrogram=True)
#sc.pl.heatmap(combined_var, genes, groupby='dataset', dendrogram=True)