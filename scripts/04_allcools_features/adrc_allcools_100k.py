import pandas as pd
import numpy as np
#import scanpy as sc
import matplotlib.pyplot as plt

from ALLCools.mcds import MCDS
#from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import *

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
var_dim = 'chrom100k'  # feature

# feature cov cutoffs
min_cov = 500
max_cov = 3000

# Regions to remove during the clustering analysis
# change this to the path to ENCODE blacklist.
# The ENCODE blacklist can be downloaded from https://github.com/Boyle-Lab/Blacklist/
black_list_path = '../../data/genome/mm10-blacklist.v2.bed.gz'
black_list_fraction = 0.2
exclude_chromosome = ['chrY','chrX','chrL']

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

mcds_path='/data/cfDNA_Jan2025/allcools/outpc.mcds_100k'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='chrom100k'
)

var_dim='chrom100k'

mcds.add_feature_cov_mean(var_dim=var_dim)

mcds.add_mc_frac(
    normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
    clip_norm_value=10000  # clip outlier values above 10 to 10
)

mcg_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type='CGN',
                                 n_top_feature=1000,
                                 plot=False)
mch_adata = mcds.get_adata(mc_type='CGN',
                           var_dim=var_dim,
                           select_hvf=True)

mch_adata.write_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom100k_features1k.h5ad')

