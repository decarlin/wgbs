import pandas as pd
import numpy as np
#import scanpy as sc
import matplotlib.pyplot as plt
import os
from ALLCools.mcds import MCDS
import gzip
from collections import defaultdict
from tqdm import tqdm
from pybedtools import BedTool

mcds_path='/data/combined_wgbs/pd_ecker.mcds'
#mcds_path='/data/combined_wgbs/ad_ecker.mcds'
#mcds_path='/data/combined_wgbs/dlb_ecker.mcds'

out_path='/data/combined_wgbs/pd_ecker100_celfie_input_10232025.txt'
#out_path='/data/combined_wgbs/ad_ecker100_celfie_input_11032025.txt'
#out_path='/data/combined_wgbs/dlb_ecker100_celfie_input_10232025.txt'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='eckerDMR'
)

var_dim='eckerDMR'

#ref_path = "/home/user/celfie/ecker_celfie_panel.tsv"
ref_path = "/home/user/celfie/ecker100_celfie_panel.tsv"

# Load reference
print("Loading reference data...")
ref = pd.read_csv(ref_path, sep='\t', header=0)

ref = ref.astype({"start": int, "end": int})

ref['chr']=ref['chr'].str.replace("^chr", "", regex=True)

#from celfie_bed.py, the names should be the same

#reconstruct the experimental data

cov = mcds['eckerDMR_da'].sel(count_type='cov')
cov_df = cov.to_pandas()

mc = mcds['eckerDMR_da'].sel(count_type='mc')
mc_df = mc.to_pandas()

# --- Extract methylated and coverage counts ---
mc_t = mc_df.T
cov_t = cov_df.T

exp_df = (
    pd.concat([mc_t, cov_t], axis=1, keys=['mc', 'cov'])
    .swaplevel(axis=1)
    .sort_index(axis=1, level=0)
)
exp_df = exp_df.reindex(['mc', 'cov'], axis=1, level=1)
exp_df.columns = [f"{cell}_{ctype}" for cell, ctype in exp_df.columns]
exp_df.columns = (
    exp_df.columns
    .str.replace("_mc$", "_METH", regex=True)
    .str.replace("_cov$", "_DEPTH", regex=True)
)
exp_df.index.name = "region"


# --- Concatenate side-by-side ---


combined_df = exp_df.join(ref, how="inner")

cols_to_prepend = ['chr', 'start', 'end']

# Duplicate the columns
dup = combined_df[cols_to_prepend].copy()

# Prepend them
combined_df = pd.concat([dup, combined_df], axis=1)

combined_df = combined_df.drop(columns='name')

# --- Write to file ---
combined_df.to_csv(out_path, sep="\t", index=False)
print("✅ Combined CELFIE input matrix written to:", out_path)

n_samples = len(mcds['cell'])
print(f"\n✅ Found {n_samples} experimental samples")
