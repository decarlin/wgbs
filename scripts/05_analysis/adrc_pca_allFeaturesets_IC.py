import scanpy as sc
import os
import pandas as pd
import numpy as np

metadata = pd.read_csv('/data/pdbp_01302024/metadata/pdbp_metadata_abbrev_0411.csv',index_col=0)

mch_adata = sc.read_h5ad(filename='/data/pdbp_01302024/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata.obs=mch_adata.obs.join(metadata)

os.chdir('/data/pdbp_01302024/allcools/')

sc.set_figure_params(dpi_save=300,format='tif')
sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_byDX',color='DX',size=50)
sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_byDX',color='DX',size=50)

sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_bySex',color='SEX',size=50)
sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_bySex',color='SEX',size=50)

mch_adata.obs['Batch..'] = mch_adata.obs['Batch..'].astype("category")

sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_byBatch',color='Batch..',size=50)
sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_byBatch',color='Batch..',size=50)


control_group = "Control"
sc.tl.rank_genes_groups(mch_adata, groupby="DX", method='t-test', batch_key="SEX", reference=control_group)


# Extract ranked gene results
result = mch_adata.uns['rank_genes_groups']
groups = [g for g in result['names'].dtype.names if g != control_group]  # Exclude control group

# Collect data
gene_data = []
for group in groups:
    for i in range(len(result['names'][group])):  # Loop through each ranked gene
        gene_name = result['names'][group][i]
        # Fetch genomic information from var
        chrom = mch_adata.var.loc[gene_name, 'chrom'] if gene_name in mch_adata.var.index else None
        start = mch_adata.var.loc[gene_name, 'start'] if gene_name in mch_adata.var.index else None
        end = mch_adata.var.loc[gene_name, 'end'] if gene_name in mch_adata.var.index else None
        gene_data.append([
            group,  # Case group
            control_group,  # Control group
            gene_name,  # Gene name
            chrom,  # Chromosome
            start,  # Start position
            end,  # End position
            result['scores'][group][i],  # Wilcoxon score
            result['logfoldchanges'][group][i],  # Log fold change
            result['pvals'][group][i],  # p-value
            result['pvals_adj'][group][i]  # Adjusted p-value
        ])

# Create DataFrame
df = pd.DataFrame(gene_data, columns=['Case_Group', 'Control_Group', 'Gene', 'Chromosome', 'Start', 'End', 'Score', 'LogFoldChange', 'PValue', 'PValueAdj'])

# Save to CSV
df.to_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_1k.csv", index=False)

#annotate with genes

import pandas as pd
import pybedtools

# Load your differentially expressed genes (DEG) data
df = pd.read_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_1k.csv")

# Convert to BED format (0-based start position for compatibility with BEDTools)
df_bed = df[['Chromosome', 'Start', 'End', 'Gene']].copy()
#df_bed['Start'] = df_bed['Start'] - 1  # Convert to 0-based BED format
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
overlaps_df.to_csv("/data/pdbp_01302024/allcools/annotated_genes.csv", index=False)
print("Annotation completed! Overlapping genes saved.")

import pandas as pd

# Load original gene list (the one used for pybedtools)
original_df = pd.read_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_1k.csv")

# Load the overlaps file produced by pybedtools
overlaps_df = pd.read_csv("/data/pdbp_01302024/allcools/annotated_genes.csv")


# Group by (Chrom, Start, End) and combine multiple genes as a semicolon-separated list
grouped_overlaps = overlaps_df.groupby(["Chromosome", "Start", "End"])["Overlap_Gene"].apply(lambda x: ";".join(x)).reset_index()

# Perform a LEFT JOIN to retain all rows from original_df
merged_df = original_df.merge(grouped_overlaps, on=["Chromosome", "Start", "End"], how="left")

# Save to CSV
merged_df.to_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_1k_annotated.csv", index=False)

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier

# Extract features (gene expression/methylation matrix) and labels (DX column)
X = mch_adata.X  # Feature matrix
DX_labels = mch_adata.obs["DX"].values  # Diagnosis labels
DX_labels_fixed = np.array([str(x) if isinstance(x, str) else "Unknown" for x in DX_labels])
DX_labels = DX_labels_fixed


# Set up 5-fold cross-validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

#del mean_auc
#del std_auc
#plt.clf()

#categories=[dx for dx in np.unique(DX_labels) if (dx != "Control" and dx != "Unknown")]

categories=["Parkinson's Disease"]

plt.figure(figsize=(8, 6))  # Initialize ROC plot

# Iterate through each category vs. "CTRL"
for category in categories:
    print(f"Processing {category} vs. Control...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "Control" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = RandomForestClassifier(n_estimators=100)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        train_accuracy = model.score(X_filtered[train_idx], y_filtered[train_idx])
        print(f"Training accuracy: {train_accuracy:.3f}")
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
plt.title("Cross-Validated ROC Curves for DX Categories vs. Control")
plt.legend()
plt.show()

output_path = "/data/pdbp_01302024/allcools/pdbp_roc_curves_100Est_PDonly.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")


#log reg

from sklearn.linear_model import LogisticRegression

del mean_auc                                                                                                                     
del std_auc                                                                                                                      
plt.clf()   

# Iterate through each category vs. "CTRL"
for category in categories:
    print(f"Processing {category} vs. Control...")
    # Binary classification: 1 = category, 0 = "CTRL"
    y = np.array([1 if dx == category else 0 if dx == "Control" else -1 for dx in DX_labels])
    mask = y != -1  
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points
    tprs, aucs = [], []
    # Perform 5-fold CV
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        model = LogisticRegression(penalty='l1',C=0.4, solver='liblinear', random_state=42)
        model.fit(X_filtered[train_idx], y_filtered[train_idx])
        train_accuracy = model.score(X_filtered[train_idx], y_filtered[train_idx])
        print(f"Training accuracy: {train_accuracy:.3f}")
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
plt.title("Cross-Validated ROC Curves for DX Categories vs. Control")
plt.legend()
plt.show()

output_path = "/data/pdbp_01302024/allcools/pdbp_roc_curves_L1LRc04_PDonly.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")


