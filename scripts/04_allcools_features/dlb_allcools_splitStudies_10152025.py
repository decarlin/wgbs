import pandas as pd
import numpy as np
#import scanpy as sc
import matplotlib.pyplot as plt
import os

from ALLCools.mcds import MCDS
#from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import *
import xarray as xr
# Basic filtering parameters. 
# These are suggesting values, cutoff maybe different for different tissue and sequencing depths.
# change this to the paths to your MCDS files, 
# ALLCools.MCDS can handle multiple MCDS files automatically

# Dimension name used to do clustering
# This corresponding to AnnData .obs and .var
obs_dim = 'cell'  # observation
var_dim = 'chrom1k'  # feature

# Regions to remove during the clustering analysis
# change this to the path to ENCODE blacklist.
# The ENCODE blacklist can be downloaded from https://github.com/Boyle-Lab/Blacklist/
black_list_path = '../../data/genome/mm10-blacklist.v2.bed.gz'
black_list_fraction = 0.2
exclude_chromosome = ['Y','X','L']

# load to memory or not
load = True

# HVF
mch_pattern = 'CHN'
mcg_pattern = 'CGN'
n_top_feature = 20000

# PC cutoff
pc_cutoff = 0.1

# KNN
knn = -1  # -1 means auto determine

# Leiden
resolution = 1

#metadata = pd.read_csv(metadata_path, index_col=0)

mcds_path='/data/combined_wgbs/dlb.mcds'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='chrom1k'
)

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

blacklist = ['ADRC3440-16','PDBP16574','PDBP16574','PDBP16676']

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

adrc_samples = metadata.loc[
    (metadata["study"] == "ADRC")
].index

adrc_samples=adrc_samples[~adrc_samples.isin(blacklist)]
adrc_samples = pd.Index(adrc_samples).astype(str)
mcds_samples = pd.Index(mcds.coords["cell"]).astype(str)
adrc_samples = adrc_samples.intersection(mcds_samples)


pdbp_samples = metadata.loc[
    (metadata["study"] == "PDBP")
].index

pdbp_samples=pdbp_samples[~pdbp_samples.isin(blacklist)]
pdbp_samples = pd.Index(pdbp_samples).astype(str)
mcds_samples = pd.Index(mcds.coords["cell"]).astype(str)
pdbp_samples = pdbp_samples.intersection(mcds_samples)

#ADRC first

var_dim='chrom1k'

mcds_adrc = mcds.sel(cell=adrc_samples)

mcds_adrc.add_feature_cov_mean(var_dim=var_dim)

import matplotlib.pyplot as plt

df = mcds_adrc.coords['chrom1k_cov_mean'].to_pandas()

os.chdir('/data/combined_wgbs/')


mcds_adrc = mcds.sel(cell=adrc_samples)

mcds_adrc.add_feature_cov_mean(var_dim=var_dim)

import matplotlib.pyplot as plt

df = mcds_adrc.coords['chrom1k_cov_mean'].to_pandas()

os.chdir('/data/combined_wgbs/')

# Each column is a sample; rows are bins across the feature
df.plot(figsize=(8, 4))
plt.xlabel('Position across feature (binned)')
plt.ylabel('Mean coverage')
plt.title('Mean CGN Coverage Across Features')
plt.tight_layout()
plt.savefig('dlb_feature_cov_mean_1kb.tiff', dpi=300)
plt.show()

old_index = mcds_adrc.coords['chrom1k'].to_pandas()

# Make a new unique index (e.g., chrom1k_0, chrom1k_1, ...)
new_index = [f"chrom1k_{i}" for i in range(len(old_index))]

# Assign this unique index back to mcds_adrc
mcds_adrc = mcds_adrc.assign_coords(chrom1k=("chrom1k", new_index))
mcds_adrc = mcds_adrc.assign_coords(original_chrom1k=("chrom1k", old_index))

mcds_adrc = mcds_adrc.filter_feature_by_cov_mean(min_cov=5)

# remove chromosomes
mcds_adrc = mcds_adrc.remove_chromosome(exclude_chromosome)


mcds_adrc.add_mc_frac(
    normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
)

#mcg_hvf = mcds_adrc.calculate_hvf_svr(var_dim=var_dim,
#                                 mc_type='CGN',
#                                 n_top_feature=1000,
#                                 plot=False)

import numpy as np

# Assuming your mcds_adrc has feature data in a variable 'cov'
cov_data = mcds_adrc['chrom1k_da_frac'].values

variances = np.var(cov_data, axis=0)

# Set a variance threshold for feature selection
variance_threshold = 0.02  # You can adjust this value

# Select features with variance above the threshold
selected_features = variances > variance_threshold

# Apply the filter
mcds_adrc_filtered = mcds_adrc.sel(chrom1k=selected_features)

mch_adata_adrc = mcds_adrc_filtered.get_adata(mc_type='CGN',
                           var_dim=var_dim)

mch_adata_adrc.write_h5ad(filename='/data/combined_wgbs/dlb_adrc_chrom1k_features1k_varfilt002.h5ad')

#PDBP                                                                         
mcds_pdbp = mcds.sel(cell=pdbp_samples)

mcds_pdbp.add_feature_cov_mean(var_dim=var_dim)

import matplotlib.pyplot as plt

df = mcds_pdbp.coords['chrom1k_cov_mean'].to_pandas()

os.chdir('/data/combined_wgbs/')

# Each column is a sample; rows are bins across the feature                                                                               
df.plot(figsize=(8, 4))
plt.xlabel('Position across feature (binned)')
plt.ylabel('Mean coverage')
plt.title('Mean CGN Coverage Across Features')
plt.tight_layout()
plt.savefig('dlb_feature_cov_mean_1kb.tiff', dpi=300)
plt.show()

old_index = mcds_pdbp.coords['chrom1k'].to_pandas()

# Make a new unique index (e.g., chrom1k_0, chrom1k_1, ...)                                                                               
new_index = [f"chrom1k_{i}" for i in range(len(old_index))]

# Assign this unique index back to mcds_pdbp                                                                                              
mcds_pdbp = mcds_pdbp.assign_coords(chrom1k=("chrom1k", new_index))
mcds_pdbp = mcds_pdbp.assign_coords(original_chrom1k=("chrom1k", old_index))

mcds_pdbp = mcds_pdbp.filter_feature_by_cov_mean(min_cov=5)

# remove chromosomes                                                                                                                      
mcds_pdbp = mcds_pdbp.remove_chromosome(exclude_chromosome)


mcds_pdbp.add_mc_frac(
    normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix                                                  
)

# Assuming your mcds_pdbp has feature data in a variable 'cov'
cov_data = mcds_pdbp['chrom1k_da_frac'].values

variances = np.var(cov_data, axis=0)

# Set a variance threshold for feature selection
variance_threshold = 0.02  # You can adjust this value

# Select features with variance above the threshold
selected_features = variances > variance_threshold

# Apply the filter
mcds_pdbp_filtered = mcds_pdbp.sel(chrom1k=selected_features)

mch_adata_pdbp = mcds_pdbp_filtered.get_adata(mc_type='CGN',
                           var_dim=var_dim)

mch_adata_pdbp.write_h5ad(filename='/data/combined_wgbs/dlb_pdbp_chrom1k_features1k_varfilt002.h5ad')
