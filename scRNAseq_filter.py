import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# MAIN DRIVER
if __name__ == "__main__":
    anndata = sc.read_h5ad(filename="anndata_with_QCMetrics.h5ad")
    sns.histplot(anndata.obs["n_genes_by_counts"])
    #sns.histplot(anndata.obs["total_counts"])
    #sns.histplot(anndata.obs["pct_counts_mt"])
    plt.show()
    #print(test)
