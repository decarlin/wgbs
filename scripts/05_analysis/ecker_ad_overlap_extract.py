import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc

metadata = pd.read_csv('~/adrc_metadata_wPatientID.csv', index_col=0)

metadata_withCancer=pd.read_csv('~/adrc_metadata_wPatientID_BCincluded.csv', index_col=0)

os.chdir('/data/cfDNA_Jan2025/allcools/')
mch_adata = sc.read_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata_wCancer= sc.read_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata_wCancer.obs.index=mch_adata_wCancer.obs.index.str.replace('_', '-')

mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata_wCancer.obs=mch_adata_wCancer.obs.join(metadata_withCancer)


# Step 1: Load the differential methylation results
#dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

with open('/home/user/AD_ecker_overlap_features.txt', 'r') as f:
    feature_list = [line.strip() for line in f]

mch_adata_selected = mch_adata_wCancer[:,feature_list]

df = pd.DataFrame(mch_adata_selected.X.toarray() if hasattr(mch_adata_selected.X, 'toarray') else mch_adata_selected.X,
                  index=mch_adata_selected.obs_names,
                  columns=mch_adata_selected.var_names)

df.to_csv('~/AD_ecker_overlap_subset.csv')

