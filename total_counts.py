from scipy.stats import poisson

def cutoff_005(adata_obj):
    mu = adata_obj.obs['total_counts'].mean()
    poisson_dist = poisson(mu=mu)

    cutoff = poisson_dist.ppf(0.99) # p-value < 0.05
    print(cutoff)
