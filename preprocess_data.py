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
    for ds in datasets:
        adatas[ds] = sc.read_10x_mtx(data_dir, prefix=ds+"_", cache=True)
    combined = ad.concat(adatas, label="dataset")
    combined.obs_names_make_unique()

    print(combined.obs['dataset'])

    # to filter cells
    sc.pp.filter_cells(combined, min_genes=200)
    sc.pp.filter_cells(combined, min_counts=1000)

    # to filter genes
    sc.pp.filter_genes(combined, min_cells=5)
    sc.pp.filter_genes(combined, min_counts=15)

    combined.var['mt'] = combined.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    print(f'Number of Genes in Combined (after filtering): {combined.n_vars}')
    print(f'Number of Cells in Combined (after filtering): {combined.n_obs}')