import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/pd_chrom1k_featuresEckerDMR_varfilt002.h5ad')

#pc1_loadings = mch_adata.obsm['X_pca'][:, 0]
#max_index = np.argmax(pc2_loadings)
#print("Max PC1 loading:", pc1_loadings[max_index])
#print("Corresponding sample:", adata.obs_names[max_index])

#pc1_series = pd.Series(pc1_loadings, index=mch_adata.obs_names)

# Sort descending to get top
#n=3
#top_pc1 = pc1_series.sort_values(ascending=False).head(n)

#blacklist with the cancer samples

blacklist = ['PDBP16676','PDBP16574','PDBP16551']
blacklist = [] 
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)

nan_indices = mch_adata.obs[pd.isna(mch_adata.obs['sex'])].index


#get the data without cancer

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_pd_EckerDMR.png',color='dx')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byDX_pd_EckerDMR.png',color='dx')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_bySex_pd_EckerDMR.png',color='sex')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_bySex_pd_EckerDMR.png',color='sex')

#sc.pl.pca_variance_ratio(mch_adata, n_pcs=50, log=True,
#                         save='pcs_1k_variancePlot_pd_bl.png')


gene_data = []
groups=['PD']


# Define the control group (adjust this to match your dataset)
control_group = "CTRL"  # Replace with the actual control label in 'DX'

# Run differential expression analysis against the control group
sc.tl.rank_genes_groups(mch_adata, groupby="dx", method='t-test', reference=control_group)

dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
dedf=dedf.set_index('names')

dedf_with_chr=dedf.join(mch_adata.var)

dedf.to_csv("/data/combined_wgbs/pd_deVsCtrl_EckerDMR_ttest.csv")

df_bed = mch_adata.var[['chrom', 'start', 'end']].copy()
df_bed['name'] = mch_adata.var.index
df_bed.to_csv("/data/combined_wgbs/pd_regions_of_interest_Ecker.bed", sep="\t", header=False, index=False)

annotation_file = "/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"

regions_bed = pybedtools.BedTool("/data/combined_wgbs/pd_regions_of_interest_Ecker.bed")
hg38_bed = pybedtools.BedTool(annotation_file)

overlaps = regions_bed.intersect(hg38_bed, wa=True, wb=True)
overlaps_df = pd.read_csv(overlaps.fn, sep="\t", header=None)
overlaps_df.columns = ['Chromosome', 'Start', 'End', 'Region', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']

overlaps_df.to_csv("/data/combined_wgbs/pd_annotated_regions_Ecker.csv", index=False)
overlaps_df=overlaps_df.set_index('Region')

df_bed_annotated=dedf_with_chr.join(overlaps_df)
df_bed_annotated=df_bed_annotated[['scores','logfoldchanges','pvals','pvals_adj','chrom','end','start','Overlap_Gene']]

df_bed_annotated.to_csv("/data/combined_wgbs/pd_deVsCtrl_ttest_Ecker_annotated.csv")

import pandas as pd
import pybedtools
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
        model = RandomForestClassifier(n_estimators=300,random_state = 1234)
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

output_path = "/data/combined_wgbs/pd_ecker_roc_curves_300Est.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")


#hyperparameter search

category='PD'

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

categories=['PD']

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

output_path = "/data/combined_wgbs/pd_roc_curves_L1c2logistic.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")

#try the same with some feature selection using the t-test

from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

del mean_auc
del std_auc
plt.clf()

K = 20

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
        X_train, X_test = X_filtered[train_idx], X_filtered[test_idx]
        y_train, y_test = y_filtered[train_idx], y_filtered[test_idx]
        # Feature selection: t-test (ANOVA F-test) on training data only
        selector = SelectKBest(score_func=f_classif, k=K)
        X_train_sel = selector.fit_transform(X_train, y_train)
        X_test_sel = selector.transform(X_test)
        model = LogisticRegression(penalty='l1',C=2, solver='liblinear', random_state=42)
        model.fit(X_train_sel, y_train)
        y_probs = model.predict_proba(X_test_sel)[:, 1]  # Probabilities for class 1
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

output_path = "/data/combined_wgbs/pd_roc_curves_L1c2logistic_20TtestFeatures.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# Step 1: Load the differential methylation results
#dmr_df = pd.read_csv('/data/cfDNA_Jan2025/allcools/DMRs_groups_vs_control_1k_annotated.csv')

dmr_df = pd.read_csv('/data/combined_wgbs/pd_deVsCtrl_ttest.csv')

# Inspect the structure of the DMRs dataframe
dmr_df.head()
dmf_df=dmr_df.set_index('names')

# Step 2: Extract relevant features (e.g., genes or CpG sites)
# Assuming the file contains a column 'gene' or 'CpG' with the features of interest

top_features = []

# Filter rows for the current group
# Sort by 'LogFoldChange' in descending order and select top 20
top_40_group = dmr_df.nsmallest(40, 'pvals').sort_values('logfoldchanges')['names']

# Select the rows from `mch_adata_wCancer` that match the features in `selected_features`
mch_adata_selected = mch_adata[:, top_40_group]

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

scaler = MinMaxScaler()
scaled_values=scaler.fit_transform(expression_df_sorted)
expression_df_sorted = pd.DataFrame(scaled_values, columns=expression_df_sorted.columns, index=expression_df_sorted.index)

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
    this_dx=mch_adata_selected.obs.loc[tick_label.get_text(),'dx']
    if pd.isna(this_dx):
        tick_label.set_color("orange")
    if this_dx=='AD':
        tick_label.set_color("blue")
    elif this_dx=='CTRL':
        tick_label.set_color("black")
    elif this_dx=='PD':
        tick_label.set_color("green")
    elif this_dx=='DLB':
        tick_label.set_color("red")    

plt.savefig('/data/combined_wgbs/pd_features_heatmap.tiff', format="tiff", dpi=300, bbox_inches="tight")

