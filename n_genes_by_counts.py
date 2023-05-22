from scipy.stats import poisson
import matplotlib.pyplot as plt
from fitter import Fitter

def cutoff_005(adata_obj):
    mu = adata_obj.obs['n_genes_by_counts'].mean()
    poisson_dist = poisson(mu=mu)

    cutoff = poisson_dist.ppf(0.9999) # p-value < 0.05

    x = range(adata_obj.obs['n_genes_by_counts'].min(), adata_obj.obs['n_genes_by_counts'].max() + 1)
    pmf = poisson_dist.pmf(x)

    plt.hist(adata_obj.obs['n_genes_by_counts'], bins=range(adata_obj.obs['n_genes_by_counts'].min(), adata_obj.obs['n_genes_by_counts'].max() + 2), 
             density=True, alpha=0.5, label='Data')
    plt.plot(x, pmf, label='Poisson Distribution')
    plt.axvline(x=cutoff, color='red', linestyle='--', label='Cutoff Value')
    plt.show()


    print(cutoff)
