import scanpy as sc
import os
import pandas as pd
import numpy as np

metadata = pd.read_csv('~/adrc_metadata_wPatientID.csv', index_col=0)

os.chdir('/data/cfDNA_Jan2025/allcools/')
mch_adata = sc.read_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom100k_features1k.h5ad')

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')

mch_adata.obs=mch_adata.obs.join(metadata)

nan_indices = mch_adata.obs[pd.isna(mch_adata.obs['DX'])].index

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_100k_PC1_PC2_byDX.png',color='DX')
sc.pl.pca(mch_adata,components='3,4',save='_100k_PC3_PC4_byDX.png',color='DX')

sc.pl.pca(mch_adata,save='_100k_PC1_PC2_bySex.png',color='SEX')
sc.pl.pca(mch_adata,components='3,4',save='_100k_PC3_PC4_bySex.png',color='SEX')

mch_adata.obs['Batch..'] = mch_adata.obs['Batch..'].astype("category")

sc.pl.pca(mch_adata,save='_100k_PC1_PC2_byBatch.png',color='Batch..')
sc.pl.pca(mch_adata,components='3,4',save='_100k_PC3_PC4_byBatch.png',color='Batch..')


sc.pl.pca_variance_ratio(mch_adata, n_pcs=50, log=True,
                         save='pcs_100k_variancePlot.png')

# Extract PCA embeddings
pc4_loadings = mch_adata.obsm["X_pca"][:, 3]  # PC4 is at index 3

# Get the indices of the two samples with the highest loadings
top_indices = np.argsort(pc4_loadings)[-2:]  # Last two are the largest

# Retrieve sample names
top_samples = mch_adata.obs.index[top_indices]

sc.tl.rank_genes_groups(mch_adata, groupby='DX', method="wilcoxon", batch_key="SEX")

df = pd.DataFrame({
    "gene": results["names"][groups[0]],  # Gene names
})

for group in groups:
    df[f"logFC_{group}"] = results["logfoldchanges"][group]
    df[f"pval_{group}"] = results["pvals"][group]
    df[f"adj_pval_{group}"] = results["pvals_adj"][group]

df.to_csv("differential_expression_full_results.csv", index=False)
