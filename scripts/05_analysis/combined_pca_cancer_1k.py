import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/chrom1k_CGN_03var.h5ad')

#pc2_loadings = mch_adata.obsm['X_pca'][:, 1]
#max_index = np.argmax(pc2_loadings)
#print("Max PC1 loading:", pc1_loadings[max_index])
#print("Corresponding sample:", adata.obs_names[max_index])

#pc2_series = pd.Series(pc2_loadings, index=mch_adata.obs_names)

# Sort descending to get top
#n=4
#top_pc2 = pc2_series.sort_values(ascending=False).head(n)

#blacklist with the cancer samples

blacklist = ['ADRC3440-16','PDBP16574','PDBP16676','PDBP16551','ADRC3533-16']



mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)

nan_indices = mch_adata.obs[pd.isna(mch_adata.obs['sex'])].index


#get the data without cancer

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_with_cancer.png',color='dx')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byDX_with_cancer.png',color='dx')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_bySex_with_cancer.png',color='sex')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_bySex_with_cancer.png',color='sex')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byStudy_with_cancer.png',color='study')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byStudy_with_cancer.png',color='study')

sc.pl.pca_variance_ratio(mch_adata, n_pcs=50, log=True,
                         save='pcs_1k_variancePlot.png')

mch_adata_cancer=mch_adata[nan_indices,]

#using combat for covariate correction
mch_adata.obs["batch_cov"] = mch_adata.obs["study"].astype(str) + "_" + mch_adata.obs["sex"].astype(str)
sc.pp.combat(mch_adata, key="batch_cov")


#DE genes

sc.tl.rank_genes_groups(mch_adata, groupby='dx', method="wilcoxon")

gene_data = []
groups=['AD','CTRL','DLB','PD']

dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
dedf.to_csv("/data/combined_wgbs/combined_de_OneVsAll.csv")

# Define the control group (adjust this to match your dataset)
control_group = "CTRL"  # Replace with the actual control label in 'DX'

# Run differential expression analysis against the control group
sc.tl.rank_genes_groups(mch_adata, groupby="dx", method="wilcoxon", reference=control_group)

dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
dedf.to_csv("/data/combined_wgbs/combined_de_OneVsCtrl.csv")


# Extract ranked gene results
result = mch_adata.uns['rank_genes_groups']
groups = [g for g in result['names'].dtype.names if g != control_group]  # Exclude control group

# Define the control group (adjust this to match your dataset)
control_group = "CTRL"  # Replace with the actual control label in 'DX'

# Run differential expression analysis against the control group
sc.tl.rank_genes_groups(mch_adata, groupby="DX", method="wilcoxon", batch_key="SEX", reference=control_group)


import pandas as pd
import pybedtools



print("✅ Merging complete! Saved as 'final_annotated_results.csv'")


import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier

# Load AnnData object (assuming it's already in memory)
# mch_adata = sc.read_h5ad("your_file.h5ad")  # Uncomment if loading from a file

# Extract features (gene expression/methylation matrix) and labels (DX column)
X = mch_adata.X  # Feature matrix
DX_labels = mch_adata.obs["dx"].values  # Diagnosis labels

# Identify unique DX categories (excluding "CTRL")
categories = [dx for dx in np.unique(DX_labels) if (dx != "CTRL" and dx != "Unknown")]

# Set up 5-fold cross-validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

#del mean_auc
#del std_auc
#plt.clf()

plt.figure(figsize=(8, 6))  # Initialize ROC plot

# Iterate through each category vs. "CTRL"
for category in categories:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = RandomForestClassifier(n_estimators=50)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR values at common FPR points
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        # Ensure the curve starts at (0,0)
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)


# Plot chance level
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/cfDNA_Jan2025/allcools/roc_curves_RF50Est.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")


#hyperparameter search

category='DLB'

y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
mask = y != -1
X_filtered, y_filtered = X[mask], y[mask]


from sklearn.model_selection import RandomizedSearchCV

rf = RandomForestClassifier(random_state=42, n_jobs=-1)

param_dist = {
    'n_estimators': [100, 300, 500, 1000],
    'max_depth': [10, 20, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 5, 10],
    'max_features': ['sqrt', 'log2']
}

random_search = RandomizedSearchCV(rf, param_dist, n_iter=20, cv=5, scoring='roc_auc', random_state=42, n_jobs=-1)
random_search.fit(X_filtered, y_filtered)

print("Best Parameters:", random_search.best_params_)

from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

param_grid = {
    'n_estimators': [100, 300, 500],
    'max_depth': [10, 20, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 5, 10],
    'max_features': ['sqrt', 'log2']
}

rf = RandomForestClassifier(random_state=42, n_jobs=-1)
grid_search = GridSearchCV(rf, param_grid, cv=5, scoring='roc_auc', verbose=1)
grid_search.fit(X_filtered, y_filtered)

print("Best Parameters:", grid_search.best_params_)

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import numpy as np
import matplotlib.pyplot as plt

del mean_auc                                                                                                                     
del std_auc                                                                                                                      
plt.clf()                                                                                                                        

# Iterate through each category vs. "CTRL"
for category in categories:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = LogisticRegression(penalty='l1',C=2, solver='liblinear', random_state=42)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR values at common FPR points
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        # Ensure the curve starts at (0,0)
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)

# Plot chance level                                                                                                              
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings                                                                                                            
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/cfDNA_Jan2025/allcools/roc_curves_L1c2logistic.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")



del mean_auc                                                                                                                     
del std_auc                                                                                                                      
plt.clf()                                                                                                                        

# Iterate through each category vs. "CTRL"
for category in ['AD']:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = LogisticRegression(penalty='l1',C=2, solver='liblinear', random_state=42)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR values at common FPR points
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        # Ensure the curve starts at (0,0)
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)

# Plot chance level                                                                                                              
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings                                                                                                            
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/cfDNA_Jan2025/allcools/roc_curves_L1c2logistic_ADonly.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")

del mean_auc                                                                                                                     
del std_auc                                                                                                                      
plt.clf()                                                                                                                        

# Iterate through each category vs. "CTRL"
for category in ['DLB']:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = LogisticRegression(penalty='l1',C=2, solver='liblinear', random_state=42)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR values at common FPR points
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        # Ensure the curve starts at (0,0)
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)

# Plot chance level                                                                                                              
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings                                                                                                            
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/cfDNA_Jan2025/allcools/roc_curves_L1c2logistic_DLBonly.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")


plt.figure(figsize=(8, 6))  # Initialize ROC plot

del mean_auc
del std_auc
plt.clf()
# Iterate through each category vs. "CTRL"
for category in ['PDD-PD']:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR values at common FPR points
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        # Ensure the curve starts at (0,0)
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)


# Plot chance level
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/cfDNA_Jan2025/allcools/roc_curves_RF100Est_PDonly.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Load the differential methylation results
dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

# Inspect the structure of the DMRs dataframe
dmr_df.head()

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the features of interest
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Load the differential methylation results
#dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/differential_expression_full_TtestResults_1k.csv')

# Inspect the structure of the DMRs dataframe
dmr_df.head()

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the features of interest

top_features = []

# Iterate over each unique 'Case_Group' or 'Clustert
for group in dmr_df['Cluster'].unique():
    # Filter rows for the current group
    group_df = dmr_df[dmr_df['Cluster'] == group]    
    # Sort by 'LogFoldChange' in descending order and select top 20
    top_20_group = group_df.nlargest(20, 'LogFoldChange')
    # Append to the list
    top_features.append(top_20_group)


# Step 4: Combine the top features for all groups into a single DataFrame
top_features_df = pd.concat(top_features)

selected_features=top_features_df['Gene']
# Select the rows from `mch_adata_wCancer` that match the features in `selected_features`
mch_adata_selected = mch_adata_wCancer[:, selected_features]

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

plt.clr()
# Step 5: Create a heatmap
plt.figure(figsize=(20, 8))
g=sns.heatmap(expression_df_sorted.T,yticklabels=False, cmap='coolwarm', annot=False, cbar_kws={'label': 'Expression Level'})
plt.title("Heatmap of Differentially Methylated Regions (DMRs) in Methylation Data")
plt.ylabel('Selected 1kb methylation windows')
plt.xlabel('Samples')
plt.show()


#for tick, label in zip(range(len(sample_annotations)), sample_annotations['DX'].fillna('Cancer')):
#    plt.text(len(expression_df_sorted.T) + 1, tick, label, color='black', fontsize=12, ha='left')

for tick_label in g.axes.get_xticklabels():
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'DX']
    if pd.isna(this_dx):
        tick_label.set_color("orange")
    if this_dx=='AD':
        tick_label.set_color("blue")
    elif this_dx=='CTRL':
        tick_label.set_color("black")
    elif this_dx=='PDD-PD':
        tick_label.set_color("green")
    elif this_dx=='DLB':
        tick_label.set_color("red")    

plt.savefig('/data/cfDNA_Jan2025/allcools/adrc_features_heatmap.tiff', format="tiff", dpi=300, bbox_inches="tight")


#celfie prep

import os
import gzip
import pandas as pd
from pybedtools import BedTool

ref_path = "/home/user/celfie/tim_matrix.txt"

# Load reference
print("Loading reference data...")
ref_cols = ['chr', 'start', 'end'] + [f"{t}_meth" for t in tissues] + [f"{t}_depth" for t in tissues]
ref = pd.read_csv(ref_path, sep='\t', header=0)
ref.columns = ref_cols
ref = ref.astype({"start": int, "end": int})

# Adjust reference to use 250 base pair windows around each CpG site
ref['start'] = ref['start'] - 250
ref['end'] = ref['end'] + 250

# Set the correct index for the reference
ref = ref.set_index(['chr', 'start', 'end'])

# Debug: Check the index of the ref DataFrame
print("Reference DataFrame index:", ref.index)

# Reset the index so 'chr', 'start', 'end' are columns
ref_reset = ref.reset_index()

# Create a BedTool from the reference intervals
bed_ref = BedTool.from_dataframe(ref_reset[['chr', 'start', 'end']])

# Build interval lookup by chromosome
window_intervals = defaultdict(list)
for chrom, start, end in zip(ref_reset['chr'], ref_reset['start'], ref_reset['end']):
    window_intervals[chrom].append((start, end))

mch_adata.var['chrom']='chr'+mch_adata.var['chrom'].astype(str)

bed_data=BedTool.from_dataframe(mch_adata.var[['chrom', 'start', 'end']])

intersected=bed_data.intersect(bed_ref)
df_inter = inter.to_dataframe()


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
bt_ref = BedTool.from_dataframe(ref.reset_index()[['chr', 'start', 'end']])

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

ref_reset_keyed = ref_reset.set_index(['chr', 'start', 'end'])
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
ref_reset_keyed.index = ref_reset_keyed.index.astype(str)

# Now perform the join
final = meth_pivot.join(depth_pivot, how='left').join(ref_reset_keyed, how='left')
final.to_csv("/data/pdbp_01302024/celfie/celfie_input.txt", sep="\t", index=True)

non_sample_cols = [col for col in final.columns if not any(x in col for x in ['_meth', '_depth'])]

# Get sample names by stripping suffixes
import re
sample_base_names = sorted(set(re.sub(r'_(meth|depth)$', '', col) 
                               for col in final.columns 
                               if '_meth' in col or '_depth' in col))

# Interleave meth and depth columns
reordered_cols = non_sample_cols + [f"{sample}_meth" for sample in sample_base_names] + [f"{sample}_depth" for sample in sample_base_names]

# Or, if you want true interleaving (meth, depth per sample side by side):
reordered_cols = non_sample_cols + [col for sample in sample_base_names for col in [f"{sample}_meth", f"{sample}_depth"]]

# Apply the order
final = final[reordered_cols]
final.columns = [col.upper() for col in final.columns]
final.to_csv("/data/pdbp_01302024/celfie/celfie_input_reordered.txt", sep="\t", index=True)



#annotate with genes

import pandas as pd
import pybedtools

# Load your differentially expressed genes (DEG) data

#df = pd.read_csv("/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k.csv")

df = pd.read_csv('/data/cfDNA_Jan2025/allcools/differential_expression_full_TtestResults_1k.csv')

# Convert to BED format (0-based start position for compatibility with BEDTools)
df_bed = df[['Chromosome', 'Start', 'End', 'Gene']].copy()
df_bed['Start'] = df_bed['Start'] - 1  # Convert to 0-based BED format
df_bed.to_csv("genes_of_interest.bed", sep="\t", header=False, index=False)

# Load hg38 annotation file (e.g., GENCODE or UCSC refGene)
annotation_file = "/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"  # Replace with the correct BED file

# Convert both files to BedTool objects
genes_bed = pybedtools.BedTool("genes_of_interest.bed")
hg38_bed = pybedtools.BedTool(annotation_file)

# Find overlapping genes
overlaps = genes_bed.intersect(hg38_bed, wa=True, wb=True)

# Convert results to DataFrame
overlaps_df = pd.read_csv(overlaps.fn, sep="\t", header=None)
overlaps_df.columns = ['Chromosome', 'Start', 'End', 'Gene', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']

# Save annotated results
overlaps_df.to_csv("/data/cfDNA_Jan2025/allcools/annotated_genes_groups.csv", index=False)
print("Annotation completed! Overlapping genes saved.")

import pandas as pd

# Load original gene list (the one used for pybedtools)
original_df = pd.read_csv("/data/cfDNA_Jan2025/allcools/differential_expression_full_TtestResults_1k.csv")

# Load the overlaps file produced by pybedtools
overlaps_df = pd.read_csv("/data/cfDNA_Jan2025/allcools/annotated_genes_groups.csv")

# Group by (Chrom, Start, End) and combine multiple genes as a semicolon-separated list
grouped_overlaps = overlaps_df.groupby(["Chromosome", "Start", "End"])["Overlap_Gene"].apply(lambda x: ";".join(x)).reset_index()


# Perform a LEFT JOIN to retain all rows from original_df
merged_df = original_df.merge(grouped_overlaps, on=["Chromosome","Start","End"], how="left")
del merged_df['column_name']

# Save to CSV
merged_df.to_csv("/data/cfDNA_Jan2025/allcools/differential_expression_full_annotated.csv", index=False)

