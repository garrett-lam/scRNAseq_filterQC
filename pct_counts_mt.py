from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns
from fitter import Fitter, get_common_distributions, get_distributions

def find_distributions(adata_obj):
    values = adata_obj.obs['pct_counts_mt'].values
    hist = sns.displot(data=values, kind="hist")
    #f = Fitter(values, distributions=['weibull_min', 'gamma', 'invgamma', 'rayleigh'])
    f = Fitter(values, distributions=get_common_distributions())
    f.fit()
    print(f.summary())
    hist.fig.savefig('bestdist_pct_counts_mt.png')
    best_dist = f.get_best()
    print()
    return best_dist