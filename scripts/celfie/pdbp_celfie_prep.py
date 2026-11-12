import scanpy as sc
import os
import pandas as pd
import numpy as np

#conda activate scanpy-env

metadata = pd.read_csv('/data/pdbp_01302024/metadata/pdbp_metadata_abbrev_0411.csv',index_col=0)

mch_adata = sc.read_h5ad(filename='/data/pdbp_01302024/allcools/chrom1k_features1k_varfilt002.h5ad')
mch_adata.obs=mch_adata.obs.join(metadata)

#celfie prep

import gzip
from pybedtools import BedTool
from collections import defaultdict

ref_path = "/home/user/celfie/tim_matrix.txt"

# Load reference
print("Loading reference data...")
ref = pd.read_csv(ref_path, sep='\t', header=0)

ref = ref.astype({"start": int, "end": int})

# Set the correct index for the reference
ref = ref.set_index(['chrom', 'start', 'end'])

# Debug: Check the index of the ref DataFrame
print("Reference DataFrame index:", ref.index)

# Reset the index so 'chr', 'start', 'end' are columns
ref_reset = ref.reset_index()

# Create a BedTool from the reference intervals
bed_ref = BedTool.from_dataframe(ref_reset[['chrom', 'start', 'end']])

# Build interval lookup by chromosome
window_intervals = defaultdict(list)
for chrom, start, end in zip(ref_reset['chrom'], ref_reset['start'], ref_reset['end']):
    window_intervals[chrom].append((start, end))

mch_adata.var['chrom']='chr'+mch_adata.var['chrom'].astype(str)

bed_data=BedTool.from_dataframe(mch_adata.var[['chrom', 'start', 'end']])

intersected=bed_data.intersect(bed_ref)
df_inter = intersected.to_dataframe()

#with values

bed_base = mch_adata.var[['chrom', 'start', 'end']].copy()
bed_base.columns = ['chrom', 'start', 'end']

# Add methylation values for all samples
meth_df = pd.DataFrame(mch_adata.X.T, columns=mch_adata.obs_names,index=mch_adata.var_names)  # shape: n_regions x n_samples
bed_meth_df = pd.concat([bed_base, meth_df], axis=1)

bed_long = bed_meth_df.melt(id_vars=['chrom', 'start', 'end'], var_name='sample', value_name='meth')

bed_long['name'] = bed_long['sample'] + ':' + bed_long['meth'].astype(str)
bed_long_clean = bed_long.dropna(subset=['chrom', 'start', 'end'])

bed_long_clean['name'] = bed_long_clean['sample'] + ':' + bed_long_clean['meth'].astype(str)
bed_long_clean['start'] = bed_long_clean['start'].astype(int)
bed_long_clean['end'] = bed_long_clean['end'].astype(int)

# Create BedTool
bt_ann = BedTool.from_dataframe(bed_long_clean[['chrom', 'start', 'end', 'name']])

# Create the reference BedTool
bt_ref = BedTool.from_dataframe(ref.reset_index()[['chrom', 'start', 'end']])

inter = bt_ann.intersect(bt_ref, wa=True, wb=True)
inter_df = inter.to_dataframe(names=['chrom', 'start', 'end', 'name', 'ref_chrom', 'ref_start', 'ref_end'])

inter_df[['sample', 'meth']] = inter_df['name'].str.split(':', expand=True)
inter_df['meth'] = inter_df['meth'].astype(float)

inter_df['depth'] = inter_df['meth'].max()

# Set index to reference regions
index_cols = ['ref_chrom', 'ref_start', 'ref_end']

# Aggregate duplicate entries by taking the mean of the methylation values
inter_df_aggregated = inter_df.groupby(index_cols + ['sample'], as_index=False)['meth'].mean()

# Pivot the DataFrame again
meth_pivot = inter_df_aggregated.pivot(index=index_cols, columns='sample', values='meth')

# Create methylation value matrix
meth_pivot = inter_df.pivot(index=index_cols, columns='sample', values='meth')

# Create dummy depth matrix

# Aggregate duplicate entries for the depth matrix by taking the mean
inter_df_depth_aggregated = inter_df.groupby(index_cols + ['sample'], as_index=False)['depth'].mean()

# Pivot the DataFrame again for the depth matrix
depth_pivot = inter_df_depth_aggregated.pivot(index=index_cols, columns='sample', values='depth')


# Rename columns for clarity
meth_pivot.columns = [f'{s}_meth' for s in meth_pivot.columns]
depth_pivot.columns = [f'{s}_depth' for s in depth_pivot.columns]


# Merge with original ref info

ref_reset_keyed = ref_reset.set_index(['chrom', 'start', 'end'])
ref_reset_keyed.index.names = ['ref_chrom', 'ref_start', 'ref_end']

meth_pivot = meth_pivot.set_index(index_cols)
depth_pivot = depth_pivot.set_index(index_cols)
ref_reset_keyed = ref_reset_keyed.set_index(index_cols)

ref_reset_keyed.index = ref_reset_keyed.index.set_levels(
    [ref_reset_keyed.index.get_level_values('ref_chrom'),
     ref_reset_keyed.index.get_level_values('ref_start').astype(int),
     ref_reset_keyed.index.get_level_values('ref_end').astype(int)],
    level=['ref_chrom', 'ref_start', 'ref_end']
)

meth_pivot.index = meth_pivot.index.astype(str)
depth_pivot.index = depth_pivot.index.astype(str)

final = meth_pivot.join(depth_pivot, how='left')

non_sample_cols = [col for col in final.columns if not any(x in col for x in ['_METH', '_DEPTH'])]

# Get sample names by stripping suffixes
import re
sample_base_names = sorted(set(re.sub(r'_(METH|DEPTH)$', '', col) 
                               for col in final.columns 
                               if '_METH' in col or '_DEPTH' in col))

# Interleave meth and depth columns
reordered_cols = non_sample_cols + [f"{sample}_meth" for sample in sample_base_names] + [f"{sample}_depth" for sample in sample_base_names]

# Or, if you want true interleaving (meth, depth per sample side by side):
reordered_cols = non_sample_cols + [col for sample in sample_base_names for col in [f"{sample}_meth", f"{sample}_depth"]]

# Apply the order
final = final[reordered_cols]
final.columns = [col.upper() for col in final.columns]
final.to_csv("/data/pdbp_01302024/celfie/pdbp_celfie_input_reordered.txt", sep="\t", index=True)
