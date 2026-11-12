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

nan_indices = mch_adata.obs[pd.isna(mch_adata.obs['DX'])].index

mch_adata = mch_adata[~mch_adata.obs['DX'].isna()].copy()

#lets do the boxplots

selected_DMR = pd.read_csv('/data/cfDNA_Jan2025/allcools/ADRC_AD_DLB_selected_DMP.txt',index_col=2,sep='\t')

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
    plt.savefig(f"/data/cfDNA_Jan2025/allcools/adrc_boxplots/{feature}_{gene_name}_boxplot_by_DX.tiff",
                dpi=300,
                format="tiff",
                bbox_inches="tight",
                pad_inches=0.1)
    plt.close()


gene_names = selected_DMR.loc[features, "Overlap_Gene"].to_list()
#write features to csv
df = pd.DataFrame(mch_adata[:, features].X.toarray(), index=mch_adata.obs_names, columns=gene_names )
df.to_csv("/data/pdbp_01302024/allcools/adrc_ad_dlb_selected_features.csv")
