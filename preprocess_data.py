import scanpy as sc
import anndata as ad
import re
import glob
import os

# preprocess directory that contains the feature, matrix, and barcodes
def preprocess_data(data_dir):
    # container for each particular dataset
    datasets = [] 
    
    # regex 
    pattern = r"(.*)_features" 

    # add all the datasets to the container
    for file in glob.glob(data_dir + '*'):
        ## get name of file without whole path
        basename = os.path.basename(file)

        ## search file and check if file is the feature file
        if re.search(pattern, basename):

            ## add name of dataset to datasets
            datasets.append(re.search(pattern,basename).group(1))

    # container to hold datasets 
    adatas = {}
    adata_obj = None

    # for each dataset read in data
    for ds in datasets:
        adatas[ds] = sc.read_10x_mtx(data_dir, prefix=ds+"_", cache=True)
    
    # combine dataset into one adata object
    adata_obj = ad.concat(adatas, label="dataset")
    adata_obj.obs_names_make_unique()

    # filter cells
    sc.pp.filter_cells(adata_obj, min_genes=200)
    sc.pp.filter_cells(adata_obj, min_counts=1000)

    # filter genes
    sc.pp.filter_genes(adata_obj, min_cells=5)
    sc.pp.filter_genes(adata_obj, min_counts=15)

    # replace mitrochrondial genes
    adata_obj.var['mt'] = adata_obj.var_names.str.startswith('MT-') 

    # calculate metrics
    sc.pp.calculate_qc_metrics(adata_obj, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # return number of cells and genes in study
    print(f'Number of Genes in adata_obj (after filtering): {adata_obj.n_vars}')
    print(f'Number of Cells in adata_obj (after filtering): {adata_obj.n_obs}')

    # return adata object
    return adata_obj