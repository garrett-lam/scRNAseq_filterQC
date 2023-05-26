from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns
from fitter import Fitter, get_common_distributions

# find ditribution function takes in the adata_obj and returns the best distribution that fits the data set
def find_distributions(adata_obj):
    # extract the n_gene_by_count values 
    values = adata_obj.obs['n_genes_by_counts'].values
    # make a histogram using the values 
    hist = sns.displot(data=values, kind="hist")
    # use the Fitter object/library to determine the best distribution that fit the data 
    f = Fitter(values, distributions=get_common_distributions())
    # plots the top 5 distributions
    f.fit()
    # prints out top 5 distributions and statistics 
    print(f.summary())
    # save histogram
    hist.fig.savefig('bestdist_n_genes_by_counts.png')
    # get the best distrubution 
    best_dist = f.get_best()
    print()
    # return best distrubution 
    return best_dist
