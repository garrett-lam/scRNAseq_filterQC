from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns
from fitter import Fitter, get_common_distributions

def cutoff_005(adata_obj):
    values = adata_obj.obs['total_counts'].values
    hist = sns.displot(data=values, kind="hist")
    f = Fitter(values, distributions=['weibull_min', 'gamma', 'invgamma', 'rayleigh'])
    f.fit()
    print(f.summary())
    hist.fig.savefig('bestdist_total_counts.png')
    print(f.get_best())