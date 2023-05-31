import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from fitter import Fitter, get_common_distributions

# find ditribution function takes in the adata_obj and returns the best distribution that fits the data set
def find_distributions(adata_obj):
     # extract the total_counts values 
    values = adata_obj.obs['total_counts'].values
    # make a histogram using the values
    hist = sns.displot(data=values, kind="hist")
    # use the Fitter object/library to determine the best distribution that fit the data 
    f = Fitter(values, distributions=get_common_distributions())
    # plots the top 5 distributions
    f.fit()
    # prints out top 5 distributions and statistics 
    print(f.summary())
    # save histogram
    hist.fig.savefig('bestdist_total_counts.png')
    # get the best distribution 
    best_dist = f.get_best()
    print()
    # return best distrubution
    return best_dist

def find_cutoff(best_dist, p_value):
    if best_dist.get('cauchy'):
        cutoff_value = stats.cauchy.ppf(q = 1-p_value, 
                                        loc=best_dist['cauchy']['loc'], 
                                        scale=best_dist['cauchy']['scale'])
        return cutoff_value
    elif best_dist.get('chi2'):
        cutoff_value = stats.chi2.ppf(q = 1-p_value, 
                                      df=best_dist['chi2']['df'], 
                                      loc=best_dist['chi2']['loc'], 
                                      scale=best_dist['chi2']['scale'])
        return cutoff_value
    elif best_dist.get('expon'):
        cutoff_value = stats.expon.ppf(q = 1-p_value, 
                                       loc=best_dist['expon']['loc'], 
                                       scale=best_dist['expon']['scale'])
        return cutoff_value
    elif best_dist.get('exponpow'):
        cutoff_value = stats.exponpow.ppf(q = 1-p_value, 
                                          b=best_dist['exponpow']['b'], 
                                          loc=best_dist['exponpow']['loc'], 
                                          scale=best_dist['exponpow']['scale'])
        return cutoff_value
    elif best_dist.get('gamma'):
        cutoff_value = stats.gamma.ppf(q = 1-p_value, 
                                       a=best_dist['gamma']['a'], 
                                       loc=best_dist['gamma']['loc'], 
                                       scale=best_dist['gamma']['scale'])
        return cutoff_value
    elif best_dist.get('lognorm'):
        cutoff_value = stats.lognorm.ppf(q = 1-p_value, 
                                         s=best_dist['lognorm']['s'], 
                                         loc=best_dist['lognorm']['loc'],
                                         scale=best_dist['lognorm']['scale'])
        return cutoff_value
    elif best_dist.get('norm'):
        cutoff_value = stats.norm.ppf(q = 1-p_value, 
                                      loc=best_dist['norm']['loc'], 
                                      scale=best_dist['norm']['scale'])
        return cutoff_value
    elif best_dist.get('powerlaw'):
        cutoff_value = stats.powerlaw.ppf(q = 1-p_value, 
                                          a=best_dist['powerlaw']['a'], 
                                          loc=best_dist['powerlaw']['loc'], 
                                          scale=best_dist['powerlaw']['scale'])
        return cutoff_value
    elif best_dist.get('rayleigh'):
        cutoff_value = stats.rayleigh.ppf(q = 1-p_value, 
                                          loc=best_dist['rayleigh']['loc'], 
                                          scale=best_dist['rayleigh']['scale'])
        return cutoff_value
    elif best_dist.get('uniform'):
        cutoff_value = stats.uniform.ppf(q = 1-p_value, 
                                         loc=best_dist['uniform']['loc'], 
                                         scale=best_dist['uniform']['scale'])
        return cutoff_value
    else:
        print('Invalid distribution')
        return 