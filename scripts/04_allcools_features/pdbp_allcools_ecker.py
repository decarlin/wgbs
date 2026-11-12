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
# To determine each cutoff more appropriately, one need to plot the distribution of each metric.
mapping_rate_cutoff = 0.5
mapping_rate_col_name = 'MappingRate'  # Name may change
final_reads_cutoff = 500000
final_reads_col_name = 'FinalReads'  # Name may change
mccc_cutoff = 0.03
mccc_col_name = 'mCCCFrac'  # Name may change
mch_cutoff = 0.2
mch_col_name = 'mCHFrac'  # Name may change
mcg_cutoff = 0.5
mcg_col_name = 'mCGFrac'  # Name may change

# change this to the paths to your MCDS files, 
# ALLCools.MCDS can handle multiple MCDS files automatically

# Dimension name used to do clustering
# This corresponding to AnnData .obs and .var
obs_dim = 'cell'  # observation
var_dim = 'eckerDMR'  # feature

# feature cov cutoffs
min_cov = 500
max_cov = 3000

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

mcds_path='/data/pdbp_01302024/allcools/pdbp_eckerDMR.mcds'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='eckerDMR'
)

var_dim='eckerDMR'

mcds.add_feature_cov_mean(var_dim=var_dim)

import matplotlib.pyplot as plt

df = mcds.coords['eckerDMR_cov_mean'].to_pandas()

os.chdir('/data/pdbp_01302024/allcools/')

# Each column is a sample; rows are bins across the feature
df.plot(figsize=(8, 4))
plt.xlabel('Position across feature (binned)')
plt.ylabel('Mean coverage')
plt.title('Mean CGN Coverage Across Features')
plt.tight_layout()
plt.savefig('pdbp_feature_cov_mean_ecker.tiff', dpi=300)
plt.show()


mcds = mcds.filter_feature_by_cov_mean(
    min_cov=5,  # minimum coverage
    max_cov=max_cov  # maximum coverage
)

# remove chromosomes
mcds = mcds.remove_chromosome(exclude_chromosome)


mcds.add_mc_frac(
    normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
)

#mcg_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
#                                 mc_type='CGN',
#                                 n_top_feature=1000,
#                                 plot=False)

import numpy as np

# Assuming your mcds has feature coverage data in a variable 'cov'
cov_data = mcds['eckerDMR_da_frac'].values

variances = np.var(cov_data, axis=0)

# Set a variance threshold for feature selection
variance_threshold = 0.01  # You can adjust this value

# Select features with variance above the threshold
selected_features = variances > variance_threshold

# Apply the filter
mcds_filtered = mcds.sel(eckerDMR=selected_features)


mch_adata = mcds_filtered.get_adata(mc_type='CGN',
                           var_dim=var_dim)

mch_adata.write_h5ad(filename='/data/pdbp_01302024/allcools/ecker_varfilt001.h5ad')

