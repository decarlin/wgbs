import scanpy as sc
import os
import pandas as pd
import numpy as np

metadata = pd.read_csv('/data/pdbp_01302024/metadata/pdbp_metadata_abbrev_0411.csv',index_col=0)

mch_adata = sc.read_h5ad(filename='/data/pdbp_01302024/allcools/chrom100k_non_cg_features1k_varfilt01.h5ad')
mch_adata.obs=mch_adata.obs.join(metadata)

os.chdir('/data/pdbp_01302024/allcools/')

#sc.set_figure_params(dpi_save=300,format='tif')
#sc.tl.pca(mch_adata)

#sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_byDX',color='DX',size=50)
#sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_byDX',color='DX',size=50)

#sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_bySex',color='SEX',size=50)
#sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_bySex',color='SEX',size=50)

mch_adata.obs['Batch..'] = mch_adata.obs['Batch..'].astype("category")
mch_adata.obs['Batch..'] = mch_adata.obs['Sex'].astype("category")

#sc.pl.pca(mch_adata,save='_pdbp_1k_PC1_PC2_byBatch',color='Batch..',size=50)
#sc.pl.pca(mch_adata,components='3,4',save='_pdbp_1k_PC3_PC4_byBatch',color='Batch..',size=50)


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
df.to_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_100k_noncg.csv", index=False)

#annotate with genes

import pandas as pd
import pybedtools

# Load your differentially expressed genes (DEG) data

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
original_df = pd.read_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_100k_noncg.csv")

# Load the overlaps file produced by pybedtools
overlaps_df = pd.read_csv("/data/pdbp_01302024/allcools/annotated_genes.csv")


# Group by (Chrom, Start, End) and combine multiple genes as a semicolon-separated list
grouped_overlaps = overlaps_df.groupby(["Chromosome", "Start", "End"])["Overlap_Gene"].apply(lambda x: ";".join(x)).reset_index()

# Perform a LEFT JOIN to retain all rows from original_df
merged_df = original_df.merge(grouped_overlaps, on=["Chromosome", "Start", "End"], how="left")

# Save to CSV
merged_df.to_csv("/data/pdbp_01302024/allcools/DMRs_groups_vs_control_100k_noncg_annotated.csv", index=False)
