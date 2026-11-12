import pandas as pd
import gzip

def load_allc_dss(
    file,
    min_cov=5,
    max_cov=100,
    cpg_only=True
):
    cols = ["chr", "pos", "strand", "context", "mc", "cov","pval"]
    
    df = pd.read_csv(
        file,
        sep="\t",
        names=cols,
        compression="gzip",
        usecols=["chr", "pos", "context", "mc", "cov"],
        dtype={
            "chr": "category",
            "pos": "int32",
            "context": "category",
            "mc": "int16",
            "cov": "int16"
        }
    )

    # --- CpG filter ---
    if cpg_only:
        df = df[df["context"].str.startswith("CG")]

    # --- Coverage filter ---
    df = df[(df["cov"] >= min_cov) & (df["cov"] <= max_cov)]

    # --- Rename for DSS ---
    df = df.rename(columns={
        "mc": "X",
        "cov": "N"
    })

    # --- Add chr prefix ---
    df["chr"] = "chr" + df["chr"].astype(str)

    return df[["chr", "pos", "N", "X"]]

df = pd.read_csv(
    "/home/user/projects/wgbs/data/sample_lists/adrc_file_list.tsv",
    sep="\t",
    header=None,
    names=["sample", "path"]
)

for _, row in df.iterrows():
    sample = row["sample"]
    path = row["path"]

    print(f"Processing {sample}...")

    d = load_allc_dss(path)

    d.to_feather(f"/data/combined_wgbs/adrc_feather/{sample}.feather")

    del d  # explicitly free memory (optional but nice)

