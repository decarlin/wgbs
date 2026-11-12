import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools
import diffxpy.api as de

metadata=pd.read_csv('/data/wgbs_pilot_2021/PilotMetadata2021.csv', index_col=2)

os.chdir('/data/combined_wgbs/pilot_plots')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/pilot_chrom1k_features1k_varfilt002.h5ad')

#pc1_loadings = mch_adata.obsm['X_pca'][:, 0]
#max_index = np.argmax(pc1_loadings)
#print("Max PC1 loading:", pc1_loadings[max_index])
#print("Corresponding sample:", mch_adata.obs_names[max_index])

#pc2_series = pd.Series(pc2_loadings, index=mch_adata.obs_names)

# Sort descending to get top
#n=3
#top_pc1 = pc1_series.sort_values(ascending=False).head(n)
#bot_pc2 = pc2_series.sort_values(ascending=True).head(n)   

#blacklist with the cancer samples

blacklist = ['KTP1','LUT1','PRT1','AD3','AD7','LC2','AD11']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]


mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_pilot_byDX.png',color='Diagnosis Group',size=30)
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_pilot_byDX.png',color='Diagnosis Group',size=30)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_pilot_bySex.png',color='Gender',size=30)
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_pilot_bySex.png',color='Gender',size=30)

