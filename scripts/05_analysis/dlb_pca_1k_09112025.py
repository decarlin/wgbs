import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools
import diffxpy.api as de

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/dlb_chrom1k_features1k_varfilt002.h5ad')

#pc1_loadings = mch_adata.obsm['X_pca'][:, 0]
#max_index = np.argmax(pc1_loadings)
#print("Max PC1 loading:", pc1_loadings[max_index])
#print("Corresponding sample:", adata.obs_names[max_index])

#pc2_series = pd.Series(pc2_loadings, index=mch_adata.obs_names)

# Sort descending to get top
#n=3
#top_pc1 = pc1_series.sort_values(ascending=False).head(n)
#bot_pc2 = pc2_series.sort_values(ascending=True).head(n)   

#blacklist with the cancer samples

blacklist = ['ADRC3440-16','PDBP16574','PDBP16574','PDBP16676']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)

#get the data without cancer

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_dlb_bl.png',color='dx')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byDX_dlb_bl.png',color='dx')

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byStudy_dlb.png',color='study')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_byStudy_dlb_bl.png',color='study')

gene_data = []
groups=['DLB']

# Define the control group (adjust this to match your dataset)
control_group = "CTRL"  # Replace with the actual control label in 'DX'

# Run differential expression analysis against the control group
#sc.tl.rank_genes_groups(mch_adata, groupby="dx", method="wilcoxon", reference=control_group)

#dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
#dedf.to_csv("/data/combined_wgbs/dlb_deVsCtrl_wilcoxon.csv")

# Run differential expression analysis against the control group

#sc.tl.rank_genes_groups(mch_adata, groupby="dx", method='logreg', covariates=['study'],reference=control_group)
sc.tl.rank_genes_groups(mch_adata, groupby="dx", method='t-test', reference=control_group)

dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
dedf=dedf.set_index('names')

dedf_with_chr=dedf.join(mch_adata.var)

dedf.to_csv("/data/combined_wgbs/dlb_deVsCtrl_ttest.csv")

df_bed = mch_adata.var[['chrom', 'start', 'end']].copy()
df_bed['name'] = mch_adata.var.index
df_bed.to_csv("/data/combined_wgbs/dlb_regions_of_interest.bed", sep="\t", header=False, index=False)

annotation_file = "/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"

regions_bed = pybedtools.BedTool("/data/combined_wgbs/dlb_regions_of_interest.bed")
hg38_bed = pybedtools.BedTool(annotation_file)

overlaps = regions_bed.intersect(hg38_bed, wa=True, wb=True)
overlaps_df = pd.read_csv(overlaps.fn, sep="\t", header=None)
overlaps_df.columns = ['Chromosome', 'Start', 'End', 'Region', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']

overlaps_df.to_csv("/data/combined_wgbs/dlb_annotated_regions.csv", index=False)
overlaps_df=overlaps_df.set_index('Region')

df_bed_annotated=dedf_with_chr.join(overlaps_df)
df_bed_annotated=df_bed_annotated[['scores','logfoldchanges','pvals','pvals_adj','chrom','end','start','Overlap_Gene']]

df_bed_annotated.to_csv("/data/combined_wgbs/dlb_deVsCtrl_ttest_annotated.csv")

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

output_path = "/data/combined_wgbs/dlb_roc_curves_300Est.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")

