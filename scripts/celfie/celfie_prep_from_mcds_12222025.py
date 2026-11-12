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

#mcds_path='/data/combined_wgbs/pd_celfie.mcds'
#mcds_path='/data/combined_wgbs/ad_celfie.mcds'
mcds_path='/data/combined_wgbs/dlb_celfie.mcds'

#out_path='/data/combined_wgbs/pd_celfie_input_10232025.txt'
#out_path='/data/combined_wgbs/ad_celfie_input_10232025.txt'
out_path='/data/combined_wgbs/dlb_celfie_input_10232025.txt'

mcds=MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='celfie'
)

var_dim='celfie'

# Constants
tissues = [
    "dendritic_cells", "endothelial_cells", "eosinophils", "erythroblasts",
    "macrophages", "monocytes", "neutrophils", "placenta", "T_cells",
    "adipose", "brain", "fibroblasts", "heart_left_ventricle", "hepatocytes",
    "lung", "mammary_gland", "megakaryocytes", "skeletal_muscle_myoblasts", "small_intestine"
]

ref_path = "/home/user/celfie/TIMs/sample_tims_summed.txt"

# Load reference
print("Loading reference data...")
ref_cols = ['chr', 'start', 'end'] + [f"{t}_meth" for t in tissues] + [f"{t}_depth" for t in tissues]
ref = pd.read_csv(ref_path, sep='\t', header=None)
ref.columns = ref_cols
ref = ref.astype({"start": int, "end": int})

ref['name'] = ref[['chr', 'start', 'end']].astype(str).agg('_'.join, axis=1)
ref['chr']=ref['chr'].str.replace("^chr", "", regex=True)

#from celfie_bed.py, the names should be the same

#reconstruct the experimental data

cov = mcds['celfie_da'].sel(count_type='cov')
cov_df = cov.to_pandas()

mc = mcds['celfie_da'].sel(count_type='mc')
mc_df = mc.to_pandas()

# --- Extract methylated and coverage counts ---
mc_t = mc_df.T
cov_t = cov_df.T

exp_df = (
    pd.concat([mc_t, cov_t], axis=1, keys=['mc', 'cov'])
    .swaplevel(axis=1)
    .sort_index(axis=1, level=0)
)
exp_df.columns = [f"{cell}_{ctype}" for cell, ctype in exp_df.columns]
exp_df.columns = (
    exp_df.columns
    .str.replace("_mc$", "_METH", regex=True)
    .str.replace("_cov$", "_DEPTH", regex=True)
)
exp_df.index.name = "region"

index_as_str = exp_df.index.astype(str)

split_parts = []
for r in index_as_str:
    parts = str(r).split("_")
    if len(parts) == 3:
        chrom, start, end = parts
    elif len(parts) == 4 and parts[0].lower().startswith("chr"):
        # handles e.g. 'chr_1_123_124'
        chrom, start, end = parts[1], parts[2], parts[3]
    else:
        raise ValueError(f"Unexpected region format: {r}")
    split_parts.append((chrom, start, end))

exp_coords = pd.DataFrame(split_parts, columns=["chrom", "start", "end"])
exp_coords["chrom"] = "chr" + exp_coords["chrom"].astype(str).str.removeprefix("chr")
exp_coords["start"] = exp_coords["start"].astype(int)
exp_coords["end"] = exp_coords["end"].astype(int)

exp_df.insert(0, "CHROM", exp_coords["chrom"])
exp_df.insert(1, "START", exp_coords["start"])
exp_df.insert(2, "END", exp_coords["end"])
exp_df = exp_df.reset_index(drop=True)


# --- Reference data ---
ref_tissue_cols = []
for t in tissues:
    ref_tissue_cols.extend([f"{t}_meth", f"{t}_depth"])

ref_df = ref[["chr", "start", "end"] + ref_tissue_cols].copy()
ref_df.columns = ["CHROM", "START", "END"] + [
    c.replace("_meth", "_METH").replace("_depth", "_DEPTH") for c in ref_tissue_cols
]

# --- Concatenate side-by-side ---
combined_df = pd.concat([exp_df, ref_df], axis=1)

# --- Write to file ---
combined_df.to_csv(out_path, sep="\t", index=False)
print("✅ Combined CELFIE input matrix written to:", out_path)

n_samples = len(mcds['cell'])
print(f"\n✅ Found {n_samples} experimental samples")
