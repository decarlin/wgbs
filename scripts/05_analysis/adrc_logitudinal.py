import scanpy as sc
import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests

metadata = pd.read_csv('~/adrc_metadata_wPatientID_BCincluded_abbrev.csv', index_col=0)

metadata_withCancer=pd.read_csv('~/adrc_metadata_wPatientID_BCincluded.csv', index_col=0)

os.chdir('/data/cfDNA_Jan2025/allcools/')
mch_adata = sc.read_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata_wCancer= sc.read_h5ad(filename='/data/cfDNA_Jan2025/allcools/chrom1k_features1k_varfilt01.h5ad')
mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata_wCancer.obs.index=mch_adata_wCancer.obs.index.str.replace('_', '-')

mch_adata.obs=mch_adata.obs.join(metadata)
mch_adata_wCancer.obs=mch_adata_wCancer.obs.join(metadata_withCancer)

results = []

# Loop over each diagnosis group
for dx in ['AD','CTRL','DLB']:
    # Subset AnnData by diagnosis
    adata_dx = mch_adata[mch_adata.obs['DX'] == dx]
    # Get all patient IDs that have both draw 1 and 2
    patients = adata_dx.obs['patient_id'].value_counts()
    paired_patients = patients[patients == 2].index
    # Create mapping: patient → draw 1 and draw 2 indices
    draw1_idx = adata_dx.obs[
        (adata_dx.obs['patient_id'].isin(paired_patients)) &
        (adata_dx.obs['DrawNumber'] == 1)
    ].index
    draw2_idx = adata_dx.obs[
        (adata_dx.obs['patient_id'].isin(paired_patients)) &
        (adata_dx.obs['DrawNumber'] == 2)
    ].index
    # Sort so patient order matches between Draw1 and Draw2
    draw1 = adata_dx[draw1_idx].to_df().groupby(adata_dx.obs.loc[draw1_idx, 'patient_id']).first().sort_index()
    draw2 = adata_dx[draw2_idx].to_df().groupby(adata_dx.obs.loc[draw2_idx, 'patient_id']).first().sort_index()
    # Perform paired t-tests feature-wise
    pvals = []
    for gene in mch_adata.var_names:
        stat, pval = ttest_rel(draw1[gene], draw2[gene])
        pvals.append(pval)
    # Store results in a DataFrame
    result_df = pd.DataFrame({
        'gene': mch_adata.var_names,
        'pval': pvals,
        'DX': dx
    })
    result_df['pval_adj'] = result_df.groupby('DX')['pval'].transform(lambda p: multipletests(p, method='fdr_bh')[1])
    logfc = (draw2.mean() - draw1.mean()).rename('logFC')
    result_df = result_df.join(logfc, on='gene')
    results.append(result_df)

# Combine all results
final_df = pd.concat(results, ignore_index=True)

genomic_info = mch_adata.var[['chrom', 'start', 'end']].copy()
genomic_info['gene'] = genomic_info.index

final_df = final_df.merge(genomic_info, on='gene')

#annotate with genes                                                                                                                                           

import pandas as pd
import pybedtools

df_bed = final_df[['chrom', 'start', 'end', 'gene']].copy()
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
overlaps_df.columns = ['chrom', 'start', 'end', 'gene', 'Overlap_Chr', 'Overlap_Start', 'Overlap_End', 'Overlap_Gene']

# Save annotated results                                                                                                                                       
overlaps_df.to_csv("/data/cfDNA_Jan2025/allcools/annotated_genes.csv", index=False)
print("Annotation completed! Overlapping genes saved.")

grouped_overlaps = overlaps_df.groupby(["chrom", "start", "end"])["Overlap_Gene"].apply(lambda x: ";".join(x)).reset_index()

grouped_overlaps['chrom']=grouped_overlaps['chrom'].astype('category')

merged_df = final_df.merge(grouped_overlaps, on=["chrom", "start", "end"], how="left")

merged_df.to_csv("/data/cfDNA_Jan2025/allcools/longitudinal_DMRs_1k_annotated.csv",index=False)
