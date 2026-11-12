import numpy as np
import scanpy as sc
import os
import pandas as pd
import scanpy.external as sce
import pybedtools
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm

# =============================================================================
# Load metadata
# =============================================================================

metadata_adrc = pd.read_csv(
    '~/adrc_metadata_wPatientID_BCincluded_abbrev.csv',
    index_col=0
)

# Match AnnData sample naming
metadata_adrc.index = 'ADRC' + metadata_adrc.index.astype(str)

metadata = pd.read_csv(
    '/data/combined_wgbs/combined_metadata.csv',
    index_col=0
)


# =============================================================================
# Create longitudinal variables
# =============================================================================

metadata_adrc["time_from_baseline"] = (
    metadata_adrc.groupby("patient_id")["AGE_DD"]
    .transform(lambda x: x - x.min())
)

metadata_adrc["baseline_age"] = (
    metadata_adrc.groupby("patient_id")["AGE_DD"]
    .transform("min")
)


# =============================================================================
# Load AnnData
# =============================================================================

os.chdir('/data/combined_wgbs/')

mch_adata = sc.read_h5ad(
    filename='/data/combined_wgbs/dlb_adrc_chrom1k_features1k_varfilt002.h5ad'
)


# =============================================================================
# Blacklist problematic samples
# =============================================================================

blacklist = [
    'ADRC3440-16',
    'ADRC3533-16'
]

mch_adata = mch_adata[
    ~mch_adata.obs_names.isin(blacklist),
    :
]


# =============================================================================
# Harmonize sample names
# =============================================================================

mch_adata.obs.index = (
    mch_adata.obs.index.str.replace('_', '-')
)

# =============================================================================
# Join metadata
# =============================================================================

mch_adata.obs = mch_adata.obs.join(metadata)
mch_adata.obs = mch_adata.obs.join(metadata_adrc)


# =============================================================================
# Sanity checks
# =============================================================================

print("Missing patient_id:",
      mch_adata.obs["patient_id"].isna().sum())

print("Missing AGE_DD:",
      mch_adata.obs["AGE_DD"].isna().sum())

print("Samples per patient:")
print(
    mch_adata.obs.groupby("patient_id")
    .size()
    .value_counts()
    .sort_index()
)


# =============================================================================
# Longitudinal mixed effects analysis
# =============================================================================

obs_df = mch_adata.obs.copy()

results = {
    "coef_time_dx": [],
    "pval_time_dx": []
}

X = mch_adata.X

interaction_term = (
    "time_from_baseline:"
    "C(dx, Treatment(reference='CTRL'))[T.DLB]"
)

print("X shape:", X.shape)
print("obs_df shape:", obs_df.shape)
print("n var names:", len(mch_adata.var_names))

X.shape[0] == obs_df.shape[0]
X.shape[1] == len(mch_adata.var_names)

keep = (
    mch_adata.obs["time_from_baseline"].notna() &
    mch_adata.obs["baseline_age"].notna() &
    mch_adata.obs["dx"].notna() &
    mch_adata.obs["sex"].notna() &
    mch_adata.obs["patient_id"].notna()
)

mch_adata = mch_adata[keep].copy()

X = mch_adata.X
obs_df = mch_adata.obs.copy()

for i, gene in enumerate(tqdm(mch_adata.var_names)):

    df = obs_df.copy()

    # Handle sparse matrix efficiently
    df["outcome"] = X[:, i].ravel()

    try:

        model = smf.mixedlm(
            "outcome ~ "
            "time_from_baseline * "
            "C(dx, Treatment(reference='CTRL')) + "
            "sex + baseline_age",

            data=df,

            groups=df["patient_id"]
        )

        result = model.fit(
            method="lbfgs",
            disp=False
        )

        coef = result.params.get(interaction_term)
        pval = result.pvalues.get(interaction_term)

    except Exception as e:

        print(f"Failed for {gene}: {e}")

        coef = np.nan
        pval = np.nan

    results["coef_time_dx"].append(coef)
    results["pval_time_dx"].append(pval)


# =============================================================================
# Compile results
# =============================================================================

results_df = pd.DataFrame(
    results,
    index=mch_adata.var_names
)


# =============================================================================
# FDR correction
# =============================================================================

_, fdrs, _, _ = multipletests(
    results_df["pval_time_dx"].fillna(1.0),
    method="fdr_bh"
)

results_df["FDR"] = fdrs


# =============================================================================
# Save primary results
# =============================================================================

results_df.to_csv(
    "/data/combined_wgbs/adrc_dlb_longitudinal_mixedlm.csv",
    sep="\t"
)


# =============================================================================
# Save region coordinates
# =============================================================================

mch_adata.var.to_csv(
    "/data/combined_wgbs/adrc_dlb_coords.csv",
    sep="\t",
    header=False,
    index=False
)


# =============================================================================
# Create BED file
# =============================================================================

df_bed = mch_adata.var[
    ['chrom', 'start', 'end']
].copy()

df_bed['name'] = mch_adata.var.index

bed_path = (
    "/data/combined_wgbs/"
    "ad_regions_of_interest.bed"
)

df_bed.to_csv(
    bed_path,
    sep="\t",
    header=False,
    index=False
)


# =============================================================================
# Annotate regions with genes
# =============================================================================

annotation_file = (
    "/data/genome/"
    "gencode.v47.basic.annotation_coding_reordered.bed"
)

regions_bed = pybedtools.BedTool(bed_path)

hg38_bed = pybedtools.BedTool(annotation_file)

overlaps = regions_bed.intersect(
    hg38_bed,
    wa=True,
    wb=True
)

overlaps_df = pd.read_csv(
    overlaps.fn,
    sep="\t",
    header=None
)

overlaps_df.columns = [
    'Chromosome',
    'Start',
    'End',
    'Region',
    'Overlap_Chr',
    'Overlap_Start',
    'Overlap_End',
    'Overlap_Gene'
]

overlaps_df.to_csv(
    "/data/combined_wgbs/ad_annotated_regions.csv",
    index=False
)

overlaps_df = overlaps_df.set_index('Region')


# =============================================================================
# Merge annotations with statistics
# =============================================================================

results_annotated = results_df.join(overlaps_df)

results_annotated = results_annotated[
    [
        'coef_time_dx',
        'pval_time_dx',
        'FDR',
        'Chromosome',
        'Start',
        'End',
        'Overlap_Gene'
    ]
]


# =============================================================================
# Save annotated results
# =============================================================================

results_annotated.to_csv(
    "/data/combined_wgbs/"
    "adrc_dlb_longitudinal_mixedlm_annotated.csv"
)

print("Done.")
