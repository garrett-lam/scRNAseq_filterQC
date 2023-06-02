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
    hist.fig.savefig('figures/best_dist_total_counts.png')
    
    # get the best distribution 
    best_dist = f.get_best()
    print()
    
    # return best distrubution
    return best_dist

def find_cutoff(best_dist, p_value):

    q = 1 - p_value
    cutoff_value = None
    
    if best_dist.get('cauchy'):
        cutoff_value = stats.cauchy.ppf(q, 
                                        loc=best_dist['cauchy']['loc'], 
                                        scale=best_dist['cauchy']['scale'])
    
    elif best_dist.get('chi2'):
        cutoff_value = stats.chi2.ppf(q, 
                                      df=best_dist['chi2']['df'], 
                                      loc=best_dist['chi2']['loc'], 
                                      scale=best_dist['chi2']['scale'])
    
    elif best_dist.get('expon'):
        cutoff_value = stats.expon.ppf(q, 
                                       loc=best_dist['expon']['loc'], 
                                       scale=best_dist['expon']['scale'])
    
    elif best_dist.get('exponpow'):
        cutoff_value = stats.exponpow.ppf(q, 
                                          b=best_dist['exponpow']['b'], 
                                          loc=best_dist['exponpow']['loc'], 
                                          scale=best_dist['exponpow']['scale'])
    
    elif best_dist.get('gamma'):
        cutoff_value = stats.gamma.ppf(q, 
                                       a=best_dist['gamma']['a'], 
                                       loc=best_dist['gamma']['loc'], 
                                       scale=best_dist['gamma']['scale'])
    
    elif best_dist.get('lognorm'):
        cutoff_value = stats.lognorm.ppf(q, 
                                         s=best_dist['lognorm']['s'], 
                                         loc=best_dist['lognorm']['loc'],
                                         scale=best_dist['lognorm']['scale'])
    
    elif best_dist.get('norm'):
        cutoff_value = stats.norm.ppf(q, 
                                      loc=best_dist['norm']['loc'], 
                                      scale=best_dist['norm']['scale'])
    
    elif best_dist.get('powerlaw'):
        cutoff_value = stats.powerlaw.ppf(q, 
                                          a=best_dist['powerlaw']['a'], 
                                          loc=best_dist['powerlaw']['loc'], 
                                          scale=best_dist['powerlaw']['scale'])
    
    elif best_dist.get('rayleigh'):
        cutoff_value = stats.rayleigh.ppf(q, 
                                          loc=best_dist['rayleigh']['loc'], 
                                          scale=best_dist['rayleigh']['scale'])
    
    elif best_dist.get('uniform'):
        cutoff_value = stats.uniform.ppf(q, 
                                         loc=best_dist['uniform']['loc'], 
                                         scale=best_dist['uniform']['scale'])
    
    else:
        print('Invalid distribution')
        return 
    
    return cutoff_value