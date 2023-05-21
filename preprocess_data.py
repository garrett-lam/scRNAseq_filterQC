import scanpy as sc
import anndata as ad
import re
import glob
import os

def preprocess_data(data_dir):
    datasets = [] 
    pattern = r"(.*)_features" 

    for file in glob.glob(data_dir + '*'):
        basename = os.path.basename(file) # gets name of file without whole path
        if re.search(pattern, basename):
            datasets.append(re.search(pattern,basename).group(1)) # adds name of dataset to datasets

    adatas = {}
    adata_obj = None

    for ds in datasets:
        adatas[ds] = sc.read_10x_mtx(data_dir, prefix=ds+"_", cache=True)
    adata_obj = ad.concat(adatas, label="dataset")
    adata_obj.obs_names_make_unique()

    print(adata_obj)

    # to filter cells
    sc.pp.filter_cells(adata_obj, min_genes=200)
    sc.pp.filter_cells(adata_obj, min_counts=1000)

    # to filter genes
    sc.pp.filter_genes(adata_obj, min_cells=5)
    sc.pp.filter_genes(adata_obj, min_counts=15)

    adata_obj.var['mt'] = adata_obj.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata_obj, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    print(f'Number of Genes in adata_obj (after filtering): {adata_obj.n_vars}')
    print(f'Number of Cells in adata_obj (after filtering): {adata_obj.n_obs}')