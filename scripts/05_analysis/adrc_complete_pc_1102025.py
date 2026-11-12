import numpy as np
import scanpy as sc
import os
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression


metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

metadata_adrc = pd.read_csv('~/adrc_metadata_wPatientID_BCincluded_abbrev.csv', index_col=0)
metadata_adrc.index = 'ADRC' + metadata_adrc.index.astype(str)

os.chdir('/data/combined_wgbs/')

mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/adrc_chrom1k_features1k_varfilt002.h5ad')

os.chdir('/data/combined_wgbs/adrc_complete_pc_1102025/')

sc.tl.pca(mch_adata)

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata.obs=mch_adata.obs.join(metadata_adrc)
mch_adata.obs['dx'] = mch_adata.obs['dx'].cat.add_categories(['Cancer'])
mch_adata.obs['dx'] = mch_adata.obs['dx'].fillna('Cancer')
mch_adata.obs['patient_id'] = mch_adata.obs['patient_id'].cat.add_categories(['DLS459P'])
mch_adata.obs['patient_id']=mch_adata.obs['patient_id'].fillna('DLS459P')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_all_adrc.png',color='dx')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byDX_all_adrc.png',color='dx')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_bySex_all_adrc.png',color='sex')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_bySex_all_adrc.png',color='sex')

sc.pl.pca_variance_ratio(mch_adata, n_pcs=50, log=True,
                         save='pcs_1k_variancePlots_all_adrc.png.png')


pc1_loadings = mch_adata.obsm['X_pca'][:, 0]
pc2_loadings = mch_adata.obsm['X_pca'][:, 1]

max_index = np.argmax(pc1_loadings)
print("Max PC1 loading:", pc1_loadings[max_index])
print("Corresponding sample:", mch_adata.obs_names[max_index])

max_index = np.argmax(pc2_loadings)
print("Max PC2 loading:", pc2_loadings[max_index])
print("Corresponding sample:", mch_adata.obs_names[max_index])

blacklist=['ADRC3440-16','ADRC3533-16']

mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]
sc.tl.pca(mch_adata)


sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_all_adrc_BLremoved.png',color='dx')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byDX_all_adrc_BLremoved.png',color='dx')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_bySex_all_adrc_BLremoved.png',color='sex')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_bySex_all_adrc_BLremoved.png',color='sex')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byPatient_all_adrc_BLremoved.png',color='patient_id')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byPatien_all_adrc_BLremoved.png',color='patient_id')

sc.pl.pca_variance_ratio(mch_adata, n_pcs=50, log=True,
                         save='pcs_1k_variancePlots_all_adrc_BLremoved.png')

min_index = np.argmin(pc1_loadings)

idx = np.argsort(arr)[:2]
print("After BL Min PC1 loading:", pc1_loadings[min_index])
print("Corresponding sample:", mch_adata.obs_names[min_index])


# Get metadata
obs = mch_adata.obs.copy()

# Get the first 10 PCs
pcs = mch_adata.obsm['X_pca'][:, :10]

# Add PCs to the dataframe for convenience
pc_df = pd.DataFrame(pcs, index=obs.index, columns=[f'PC{i+1}' for i in range(10)])
obs = pd.concat([obs, pc_df], axis=1)

# Keep only patients with exactly two samples
patient_counts = obs['patient_id'].value_counts()
patients_two_samples = patient_counts[patient_counts == 2].index
obs_two = obs[obs['patient_id'].isin(patients_two_samples)]

# Function to compute Euclidean distance between two rows
def euclidean_distance(row1, row2):
    return np.linalg.norm(row1 - row2)

# Calculate distances
distances = {}
for patient in patients_two_samples:
    patient_rows = obs_two.loc[obs_two['patient_id'] == patient, [f'PC{i+1}' for i in range(10)]].values
    distances[patient] = euclidean_distance(patient_rows[0], patient_rows[1])

# Convert to a Series if you want
distances_series = pd.Series(distances)

# Create a DataFrame from the distances
dist_df = pd.DataFrame({
    'patient_id': list(distances.keys()),
    'euclidean_distance': list(distances.values())
})

# Add the DX annotation for each patient
# We'll just take the first sample's DX since both samples should match
dx_map = obs_two.groupby('patient_id')['DX'].first()
dist_df['DX'] = dist_df['patient_id'].map(dx_map)

# Save to CSV
dist_df.to_csv('/data/combined_wgbs/adrc_complete_pc_1102025/patient_pairwise_distances.csv', index=False)
