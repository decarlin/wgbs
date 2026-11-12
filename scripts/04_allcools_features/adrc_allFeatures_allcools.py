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
var_dim = 'chrom1k'  # feature
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

mcds_path='/data/ADRC_fastq/allcools/outpc.mcds'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='chrom1k'
)

var_dim='chrom1k'

mcds.add_feature_cov_mean(var_dim=var_dim)

import matplotlib.pyplot as plt

df = mcds.chrom1k_cov_mean.to_pandas()

os.chdir('/data/ADRC_fastq/allcools/')

# Each column is a sample; rows are bins across the feature
df.plot(figsize=(8, 4))
plt.xlabel('Position across feature (binned)')
plt.ylabel('Mean coverage')
plt.title('Mean CGN Coverage Across Features')
plt.tight_layout()
plt.savefig('adrc_feature_cov_mean_1kb.tiff', dpi=300)
plt.show()


mcds = mcds.filter_feature_by_cov_mean(
    min_cov=5,  # minimum coverage
    max_cov=10000  # maximum coverage
)


old_index = mcds.coords['chrom1k'].to_pandas()

# Make a new unique index (e.g., chrom1k_0, chrom1k_1, ...)
new_index = [f"chrom1k_{i}" for i in range(len(old_index))]

# Assign this unique index back to mcds
mcds = mcds.assign_coords(chrom1k=("chrom1k", new_index))
mcds = mcds.assign_coords(original_chrom1k=("chrom1k", old_index))

mcds = mcds.filter_feature_by_cov_mean(min_cov=5)

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
cov_data = mcds['chrom1k_da_frac'].values

variances = np.var(cov_data, axis=0)

# Set a variance threshold for feature selection
variance_threshold = 0.02  # You can adjust this value

# Select features with variance above the threshold
selected_features = variances > variance_threshold

# Apply the filter
mcds_filtered = mcds.sel(chrom1k=selected_features)


mch_adata = mcds_filtered.get_adata(mc_type='CGN',
                           var_dim=var_dim)

mch_adata.write_h5ad(filename='/data/ADRC_fastq/allcools/chrom1k_features1k_varfilt002.h5ad')

