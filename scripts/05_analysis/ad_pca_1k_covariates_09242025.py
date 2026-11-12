import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools
from sklearn.model_selection import StratifiedGroupKFold
import statsmodels.formula.api as smf
import statsmodels.api as sm
from tqdm import tqdm  # nice progress bar
from statsmodels.stats.multitest import multipletests


metadata_adrc = pd.read_csv('~/adrc_metadata_wPatientID_BCincluded_abbrev.csv', index_col=0)
metadata_adrc.index = 'ADRC' + metadata_adrc.index.astype(str)

metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/ad_chrom1k_features1k_varfilt002.h5ad')

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


#try this with all data but with the GEE

obs_df = mch_adata.obs.copy()

gee_results = {
    "coef_dx": [],
    "pval_dx": []
}

# Loop through features
for gene in tqdm(mch_adata.var_names):
    df = obs_df.copy()
    df["outcome"] = (
        mch_adata[:, gene].X.toarray().flatten()
    )
    gee_model = smf.gee(
        "outcome ~ C(dx,Treatment(reference='CTRL')) + sex",
        groups="patient_id",
        data=df,
        family=sm.families.Gaussian(),
        cov_struct=sm.cov_struct.Exchangeable()
        )
    result = gee_model.fit()
    gee_results["coef_dx"].append(result.params.get("C(dx, Treatment(reference='CTRL'))[T.AD]"))
    gee_results["pval_dx"].append(result.pvalues.get("C(dx, Treatment(reference='CTRL'))[T.AD]"))

# Convert to DataFrame for processing
results_df = pd.DataFrame(gee_results, index=mch_adata.var_names)

# Add logFC column (same as coef_intervention, assuming Gaussian link)
#T.CTRL points the wrong way, flip it
#results_df["logFC_dx"] = -results_df["coef_dx"]

_, fdrs, _, _ = multipletests(results_df["pval_dx"].fillna(1.0), method="fdr_bh")
results_df["FDR_intervention"] = fdrs

results_df.to_csv("/data/combined_wgbs/ad_deVsCtrl_gee.csv", sep="\t", header=False, index=False)

mch_adata.var.to_csv("/data/combined_wgbs/ad_coords.csv", sep="\t", header=False, index=False)

df_bed = mch_adata.var[['chrom', 'start', 'end']].copy()
df_bed['name'] = mch_adata.var.index
df_bed.to_csv("/data/combined_wgbs/ad_regions_of_interest.bed", sep="\t", header=False, index=False)

annotation_file = "/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"

regions_bed = pybedtools.BedTool("/data/combined_wgbs/ad_regions_of_interest.bed")
hg38_bed = pybedtools.BedTool(annotation_file)

overlaps = regions_bed.intersect(hg38_bed, wa=True, wb=True)
overlaps_df = pd.read_csv(overlaps.fn, sep="\t", header=None)
overlaps_df.columns = ['Chromosome', 'Start', 'End', 'Region', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']
overlaps_df.to_csv("/data/combined_wgbs/ad_annotated_regions.csv", index=False)
overlaps_df=overlaps_df.set_index('Region')


df_bed_annotated=results_df.join(overlaps_df)
df_bed_annotated=df_bed_annotated[['coef_dx','pval_dx','FDR_intervention','Chromosome','Start','End','Overlap_Gene']]

df_bed_annotated.to_csv("/data/combined_wgbs/ad_deVsCtrl_gee_annotated.csv")
