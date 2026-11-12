import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools

metadata_adrc = pd.read_csv('~/adrc_metadata_wPatientID_BCincluded_abbrev.csv', index_col=0)
metadata_adrc.index = 'ADRC' + metadata_adrc.index.astype(str)

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/adrc_eckerDMR_varfilt002.h5ad')

#pc1_loadings = mch_adata.obsm['X_pca'][:, 0]
#pc2_loadings = mch_adata.obsm['X_pca'][:, 1]
#max_index = np.argmax(pc2_loadings)
#print("Max PC1 loading:", pc1_loadings[max_index])
#print("Corresponding sample:", adata.obs_names[max_index])

#pc1_series = pd.Series(pc1_loadings, index=mch_adata.obs_names)
#pc2_series = pd.Series(pc2_loadings, index=mch_adata.obs_names)

# Sort descending to get top
#n=1
#top_pc1 = pc1_series.sort_values(ascending=False).head(n)
#top_pc2 = pc2_series.sort_values(ascending=False).head(1)

#blacklist with the cancer samples
#


blacklist = ['ADRC3440-16','ADRC3533-16']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata.obs=mch_adata.obs.join(metadata_adrc)


#loop through different regions this with all data but with a subset of the data, pulled from the feature lists given by

obs_df = mch_adata.obs.copy()

import os
import glob
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind


region_dir = "/home/user/projects/wgbs/data/reference/ecker_regions_tables"
region_tables = glob.glob(os.path.join(region_dir, "*.tsv"))

output_dir = "/data/combined_wgbs/pca_region_plots"
os.makedirs(output_dir, exist_ok=True)

for table_path in region_tables:
    
    region_name = os.path.basename(table_path).replace("_crit01_100.tsv", "")
    print(f"\nRunning PCA for {region_name}")
    
    # read region table
    df = pd.read_csv(table_path, sep="\t")
    features = df["name"].astype(str)
    
    # intersect features
    valid_features = mch_adata.var_names.intersection(features)
    print(f"{region_name}: {len(valid_features)} valid features")
    
    if len(valid_features) < 5:
        print("Skipping — too few features")
        continue
    
    # subset features
    sub = mch_adata[:, valid_features].copy()
    
    # keep only AD and CTRL
    sub = sub[sub.obs["dx"].isin(["AD", "CTRL"]), :].copy()
    
    # scale + PCA
    sc.pp.scale(sub)
    sc.tl.pca(sub, n_comps=5)
    
    # extract PCs
    pcs = pd.DataFrame(
        sub.obsm["X_pca"][:, 0:2],
        index=sub.obs_names,
        columns=["PC1", "PC2"]
    )
    
    pcs["dx"] = sub.obs["dx"]
    
    # variance explained
    var_ratio = sub.uns["pca"]["variance_ratio"]
    pc1_var = round(var_ratio[0] * 100, 2)
    pc2_var = round(var_ratio[1] * 100, 2)
    
    # plot
    plt.figure(figsize=(6, 5))
    
    sns.scatterplot(
        data=pcs,
        x="PC1",
        y="PC2",
        hue="dx",
        palette={"AD": "firebrick", "CTRL": "steelblue"},
        s=80
    )
    
    plt.title(f"{region_name} PCA")
    plt.xlabel(f"PC1 ({pc1_var}% var)")
    plt.ylabel(f"PC2 ({pc2_var}% var)")
    plt.legend(title="Diagnosis")
    plt.tight_layout()
    
    # save
    save_path = os.path.join(output_dir, f"{region_name}_PCA.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    
    print(f"Saved: {save_path}")

    pc1_ad = pcs.loc[pcs["dx"] == "AD", "PC1"]
    pc1_ctrl = pcs.loc[pcs["dx"] == "CTRL", "PC1"]

    t_stat, p_val = ttest_ind(pc1_ad, pc1_ctrl)
    print(f"{region_name} PC1 AD vs CTRL p-value: {p_val}")
    
print("\nAll regions complete.")


