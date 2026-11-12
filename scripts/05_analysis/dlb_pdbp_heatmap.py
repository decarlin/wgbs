import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata =sc.read_h5ad(filename='/data/combined_wgbs/dlb_chrom1k_features1k_varfilt002.h5ad')

mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata=mch_adata[mch_adata.obs['study']=='PDBP',:]

blacklist = ['ADRC3440-16','PDBP16574','PDBP16574','PDBP16676']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

# Step 1: Load the differential methylation results dmr_df =
#pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

# Inspect the structure of the DMRs dataframe

# Step 1: Load the differential methylation results dmr_df =
#pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

dmr_df =
pd.read_csv('/data/combined_wgbs/dlb_pdbp_deVsCtrl_ttest_annotated.csv')

# Inspect the structure of the DMRs dataframe dmr_df.head()

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the
# features of interest

top_features = []

# Iterate over each unique 'Case_Group' or 'Clustert swithto just dlb

# Sort by 'LogFoldChange' in descending order and select top 20
bottom_20_group = dmr_df.nsmallest(20, 'logfoldchanges')
top_20_group= dmr_df.nlargest(20, 'logfoldchanges') # Append to the list
top_features.append(top_20_group) top_features.append(bottom_20_group)


# Step 4: Combine the top features for all groups into a single
DataFrame top_features_df = pd.concat(top_features)

selected_features=top_features_df['names'] # Select the rows from
`mch_adata_wCancer` that match the features in `selected_features`
mch_adata_selected = mch_adata[:, selected_features]

# Step 4: Extract the expression matrix expression_data =
mch_adata_selected.X # This gets the expression matrix (likely sparse)

# Convert to dense array if it's sparse if isinstance(expression_data,
np.ndarray) == False: expression_data = expression_data.toarray()

# Create a DataFrame for better visualization expression_df =
pd.DataFrame(expression_data, columns=mch_adata_selected.var_names,
index=mch_adata_selected.obs_names)

# Step 5: Create a heatmap sample_order =
mch_adata_selected.obs['dx'].sort_values().index sample_annotations =
mch_adata_selected.obs[['dx']].loc[sample_order] expression_df_sorted
= expression_df.loc[sample_order]

# Step 5: Create a heatmap plt.figure(figsize=(20, 8))
g=sns.heatmap(expression_df_sorted.T,yticklabels=False,
cmap='coolwarm', annot=False, cbar_kws={'label': 'Expression Level'})
plt.title("Heatmap of PDBP DLB Differentially methylated 1k regions")
plt.ylabel('Selected 1kb methylation windows') plt.xlabel('Samples')
plt.show()

#for tick, label in zip(range(len(sample_annotations)),
#    sample_annotations['DX'].fillna('Cancer')):
#    plt.text(len(expression_df_sorted.T) + 1, tick, label,
#    color='black', fontsize=12, ha='left')

for tick_label in g.axes.get_xticklabels():
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'dx'] if
    pd.isna(this_dx): tick_label.set_color("orange") elif
    this_dx=='CTRL': tick_label.set_color("black") elif
    this_dx=="Parkinson's Disease": tick_label.set_color("green") elif
    this_dx=='DLB': tick_label.set_color("red")

plt.savefig('/data/pdbp_01302024/allcools/dlb_pdbp_dmp_features_heatmap.tiff',
format="tiff", dpi=300, bbox_inches="tight")

#ADRC

os.chdir('/data/combined_wgbs/')
mch_adata =sc.read_h5ad(filename='/data/combined_wgbs/dlb_chrom1k_features1k_varfilt002.h5ad')


mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata=mch_adata[mch_adata.obs['study']=='ADRC',:]

blacklist = ['ADRC3440-16','PDBP16574','PDBP16574','PDBP16676']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

# Step 1: Load the differential methylation results dmr_df =
#pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

# Inspect the structure of the DMRs dataframe

# Step 1: Load the differential methylation results dmr_df =
#pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

dmr_df = pd.read_csv('/data/combined_wgbs/dlb_adrc_deVsCtrl_ttest_annotated.csv')

# Inspect the structure of the DMRs dataframe dmr_df.head()

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the
# features of interest

top_features = []

# Iterate over each unique 'Case_Group' or 'Clustert swithto just dlb

# Sort by 'LogFoldChange' in descending order and select top 20
bottom_20_group = dmr_df.nsmallest(20, 'logfoldchanges')
top_20_group= dmr_df.nlargest(20, 'logfoldchanges') # Append to the list
top_features.append(top_20_group)
top_features.append(bottom_20_group)


# Step 4: Combine the top features for all groups into a single DataFrame                                                                        
top_features_df = pd.concat(top_features)

selected_features=top_features_df['names']
# Select the rows from `mch_adata_wCancer` that match the features in `selected_features`                                                        
mch_adata_selected = mch_adata[:, selected_features]

# Step 4: Extract the expression matrix                                                                                                          
expression_data = mch_adata_selected.X  # This gets the expression matrix (likely sparse)                                                        

# Convert to dense array if it's sparse                                                                                                          
if isinstance(expression_data, np.ndarray) == False:
    expression_data = expression_data.toarray()

# Create a DataFrame for better visualization                                                                                                    
expression_df = pd.DataFrame(expression_data, columns=mch_adata_selected.var_names, index=mch_adata_selected.obs_names)

# Step 5: Create a heatmap                                                                                                                       
sample_order = mch_adata_selected.obs['dx'].sort_values().index
sample_annotations = mch_adata_selected.obs[['dx']].loc[sample_order]
expression_df_sorted = expression_df.loc[sample_order]

# Step 5: Create a heatmap                                                                                                                       
plt.figure(figsize=(20, 8))
g=sns.heatmap(expression_df_sorted.T,yticklabels=False, cmap='coolwarm', annot=False, cbar_kws={'label': 'Expression Level'})
plt.title("Heatmap of ADRC DLB Differentially methylated 1k regions")
plt.ylabel('Selected 1kb methylation windows')
plt.xlabel('Samples')
plt.show()

#for tick, label in zip(range(len(sample_annotations)), sample_annotations['DX'].fillna('Cancer')):                                              
#    plt.text(len(expression_df_sorted.T) + 1, tick, label, color='black', fontsize=12, ha='left')                                               

for tick_label in g.axes.get_xticklabels():
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'dx']
    if pd.isna(this_dx):
        tick_label.set_color("orange")
    elif this_dx=='CTRL':
        tick_label.set_color("black")
    elif this_dx=="Parkinson's Disease":
        tick_label.set_color("green")
    elif this_dx=='DLB':
        tick_label.set_color("red")

plt.savefig('/data/pdbp_01302024/allcools/dlb_adrc_dmp_features_heatmap.tiff', format="tiff", dpi=300, bbox_inches="tight")



