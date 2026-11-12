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
mch_adata=mch_adata[mch_adata.obs['dx']=='CTRL',:]

#get the data without cancer

sc.tl.pca(mch_adata)

sc.pl.pca(mch_adata,save='_1k_PC1_PC2_byDX_controlBothStudy.png',color='study')
sc.pl.pca(mch_adata,components='3,4',save='_1k_PC3_PC4_controlBothStudy.png',color='study')

gene_data = []
groups=['PDBP']

# Define the control group (adjust this to match your dataset)
control_group = "ADRC"  # Replace with the actual control label in 'DX'

# Run differential expression analysis against the control group
#sc.tl.rank_genes_groups(mch_adata, groupby="dx", method="wilcoxon", reference=control_group)

#dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
#dedf.to_csv("/data/combined_wgbs/dlb_pdbp_deVsCtrl_wilcoxon.csv")

# Run differential expression analysis against the control group

#sc.tl.rank_genes_groups(mch_adata, groupby="dx", method='logreg', covariates=['study'],reference=control_group)
sc.tl.rank_genes_groups(mch_adata, groupby="study", method='t-test', reference=control_group)

dedf = sc.get.rank_genes_groups_df(mch_adata, group=None)
dedf=dedf.set_index('names')

dedf_with_chr=dedf.join(mch_adata.var)

dedf.to_csv("/data/combined_wgbs/control_adrcVsPdbp_ttest.csv")

df_bed = mch_adata.var[['chrom', 'start', 'end']].copy()
df_bed['name'] = mch_adata.var.index
df_bed.to_csv("/data/combined_wgbs/dlb_pdbp_regions_of_interest.bed", sep="\t", header=False, index=False)

annotation_file = "/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"

regions_bed = pybedtools.BedTool("/data/combined_wgbs/dlb_pdbp_regions_of_interest.bed")
hg38_bed = pybedtools.BedTool(annotation_file)

overlaps = regions_bed.intersect(hg38_bed, wa=True, wb=True)
overlaps_df = pd.read_csv(overlaps.fn, sep="\t", header=None)
overlaps_df.columns = ['Chromosome', 'Start', 'End', 'Region', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']

overlaps_df.to_csv("/data/combined_wgbs/dlb_pdbp_annotated_regions.csv", index=False)
overlaps_df=overlaps_df.set_index('Region')

df_bed_annotated=dedf_with_chr.join(overlaps_df)
df_bed_annotated=df_bed_annotated[['scores','logfoldchanges','pvals','pvals_adj','chrom','end','start','Overlap_Gene']]

df_bed_annotated.to_csv("/data/combined_wgbs/control_adrcVsPdbp_ttest_annotated.csv")
