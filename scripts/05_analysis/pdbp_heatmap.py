import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc

metadata = pd.read_csv('/data/pdbp_01302024/metadata/pdbp_metadata_abbrev_0411.csv',index_col=0)

mch_adata = sc.read_h5ad(filename='/data/pdbp_01302024/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata.obs=mch_adata.obs.join(metadata)

os.chdir('/data/pdbp_01302024/allcools/')

# Step 1: Load the differential methylation results
#dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

# Inspect the structure of the DMRs dataframe

# Step 1: Load the differential methylation results
#dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

dmr_df = pd.read_csv('/data/pdbp_01302024/allcools/DMRs_groups_vs_control_1k_annotated.csv')

# Inspect the structure of the DMRs dataframe
dmr_df.head()

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the features of interest

top_features = []

# Iterate over each unique 'Case_Group' or 'Clustert
for group in dmr_df['Case_Group'].unique():
    # Filter rows for the current group
    group_df = dmr_df[dmr_df['Case_Group'] == group]    
    # Sort by 'LogFoldChange' in descending order and select top 20
    top_20_group = group_df.nsmallest(20, 'LogFoldChange')
    # Append to the list
    top_features.append(top_20_group)


# Step 4: Combine the top features for all groups into a single DataFrame
top_features_df = pd.concat(top_features)

selected_features=top_features_df['Gene']
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
sample_order = mch_adata_selected.obs['DX'].sort_values().index
sample_annotations = mch_adata_selected.obs[['DX']].loc[sample_order]
expression_df_sorted = expression_df.loc[sample_order]

# Step 5: Create a heatmap
plt.figure(figsize=(20, 8))
g=sns.heatmap(expression_df_sorted.T,yticklabels=False, cmap='coolwarm', annot=False, cbar_kws={'label': 'Expression Level'})
plt.title("Heatmap of PDBP Differentially hypomethylated Regions (DMRs) in Methylation Data")
plt.ylabel('Selected 1kb methylation windows')
plt.xlabel('Samples')
plt.show()

#for tick, label in zip(range(len(sample_annotations)), sample_annotations['DX'].fillna('Cancer')):
#    plt.text(len(expression_df_sorted.T) + 1, tick, label, color='black', fontsize=12, ha='left')

for tick_label in g.axes.get_xticklabels():
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'DX']
    if pd.isna(this_dx):
        tick_label.set_color("orange")
    elif this_dx=='Control':
        tick_label.set_color("black")
    elif this_dx=="Parkinson's Disease":
        tick_label.set_color("green")
    elif this_dx=='Dementia with Lewy Bodies':
        tick_label.set_color("red")    

plt.savefig('/data/pdbp_01302024/allcools/pdbp_hypo_features_heatmap.tiff', format="tiff", dpi=300, bbox_inches="tight")

del mch_adata_selected
del expression_data
del expression_df
del sample_order
del sample_annotations
del expression_df_sorted
del top_features

top_features = []

# Iterate over each unique 'Case_Group' or 'Clustert
for group in dmr_df['Case_Group'].unique():
    # Filter rows for the current group
    group_df = dmr_df[dmr_df['Case_Group'] == group]    
    # Sort by 'LogFoldChange' in descending order and select top 20
    top_20_group = group_df.nlargest(20, 'LogFoldChange')
    # Append to the list
    top_features.append(top_20_group)


# Step 4: Combine the top features for all groups into a single DataFrame
top_features_df = pd.concat(top_features)

selected_features=top_features_df['Gene']
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
sample_order = mch_adata_selected.obs['DX'].sort_values().index
sample_annotations = mch_adata_selected.obs[['DX']].loc[sample_order]
expression_df_sorted = expression_df.loc[sample_order]

plt.figure().clear()
# Step 5: Create a heatmap
plt.figure(figsize=(20, 8))
g=sns.heatmap(expression_df_sorted.T,yticklabels=False, cmap='coolwarm', annot=False, cbar_kws={'label': 'Expression Level'})
plt.title("Heatmap of PDBP Differentially hypermethylated Regions (DMRs) in Methylation Data")
plt.ylabel('Selected 1kb methylation windows')
plt.xlabel('Samples')
plt.show()

#for tick, label in zip(range(len(sample_annotations)), sample_annotations['DX'].fillna('Cancer')):
#    plt.text(len(expression_df_sorted.T) + 1, tick, label, color='black', fontsize=12, ha='left')

for tick_label in g.axes.get_xticklabels():
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'DX']
    if pd.isna(this_dx):
        tick_label.set_color("orange")
    elif this_dx=='Control':
        tick_label.set_color("black")
    elif this_dx=="Parkinson's Disease":
        tick_label.set_color("green")
    elif this_dx=='Dementia with Lewy Bodies':
        tick_label.set_color("red")    

plt.savefig('/data/pdbp_01302024/allcools/pdbp_hyper_features_heatmap.tiff', format="tiff", dpi=300, bbox_inches="tight")

#lets do the boxplots

selected_DMR = pd.read_csv('/data/pdbp_01302024/allcools/PDBP_PD_selected_DMP.txt',index_col=2,sep='\t')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Example: list of features you care about (index of your dataframe)
features = selected_DMR.index.tolist()

# Loop through each feature and make a boxplot
for feature in features:
    # Extract expression for the feature
    gene_name = selected_DMR.loc[feature, "Overlap_Gene"]
    data = pd.DataFrame({
        "methylation": mch_adata[:, feature].X.toarray().flatten() if hasattr(mch_adata[:, feature].X, 'toarray') else mch_adata[:, feature].X.flatten(),
        "DX": mch_adata.obs["DX"].values
    })
    plt.figure(figsize=(6, 4))
    sns.boxplot(x="DX", y="methylation", data=data)
    sns.stripplot(x="DX", y="methylation", data=data, color="black", alpha=0.3, jitter=True)
    plt.title(f"Boxplot of {feature} ({gene_name}) by DX")
    plt.tight_layout()
    plt.figure(figsize=(5, 5))
    plt.savefig(f"/data/pdbp_01302024/allcools/pdbp_boxplots/{feature}_{gene_name}_boxplot_by_DX_DMR_PD.png",
                dpi=300,
                format="tiff",
                bbox_inches="tight",
                pad_inches=0.1)
    plt.close()


gene_names = selected_DMR.loc[features, "Overlap_Gene"].to_list()
#write features to csv
df = pd.DataFrame(mch_adata[:, features].X.toarray(), index=mch_adata.obs_names, columns=gene_names )
df.to_csv("/data/pdbp_01302024/allcools/pdbp_pd_selected_features.csv")
