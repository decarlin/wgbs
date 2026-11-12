import os
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from tqdm import tqdm
from pybedtools import BedTool

# Constants
tissues = [
    "dendritic_cells", "endothelial_cells", "eosinophils", "erythroblasts",
    "macrophages", "monocytes", "neutrophils", "placenta", "T_cells",
    "adipose", "brain", "fibroblasts", "heart_left_ventricle", "hepatocytes",
    "lung", "mammary_gland", "megakaryocytes", "skeletal_muscle_myoblasts", "small_intestine"
]

ref_path = "/home/user/celfie/TIMs/sample_tims_summed.txt"
allc_dir = "/data/pdbp_01302024/aligned_hisat/allc_files_TEST"
out_path = "/data/pdbp_01302024/aligned_hisat/celfie_input_matrix_TEST.tsv"

# Load reference
print("Loading reference data...")
ref_cols = ['chr', 'start', 'end'] + [f"{t}_meth" for t in tissues] + [f"{t}_depth" for t in tissues]
ref = pd.read_csv(ref_path, sep='\t', header=None)
ref.columns = ref_cols
ref = ref.astype({"start": int, "end": int})

# Adjust reference to use 250 base pair windows around each CpG site
ref['start'] = ref['start'] - 250
ref['end'] = ref['end'] + 250

# Set the correct index for the reference
ref = ref.set_index(['chr', 'start', 'end'])

# Debug: Check the index of the ref DataFrame
print("Reference DataFrame index:", ref.index)

# Reset the index so 'chr', 'start', 'end' are columns
ref_reset = ref.reset_index()

# Create a BedTool from the reference intervals
bed_ref = BedTool.from_dataframe(ref_reset[['chr', 'start', 'end']])

# Build interval lookup by chromosome
window_intervals = defaultdict(list)
for chrom, start, end in zip(ref_reset['chr'], ref_reset['start'], ref_reset['end']):
    window_intervals[chrom].append((start, end))

# Function to stream and process a single ALLC file
def process_allc_file(filepath):
    sample_id = os.path.basename(filepath).split("_")[1].split(".")[0]
    window_meth = defaultdict(int)
    window_total = defaultdict(int)

    try:
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                chrom, pos, strand, mc, total, context = parts[:6]

                try:
                    chrom = "chr" + chrom if not chrom.startswith("chr") else chrom
                    pos = int(pos)
                    mc = int(mc)
                    total = int(total)
                except ValueError:
                    continue

                if chrom not in window_intervals:
                    continue

                for start, end in window_intervals[chrom]:
                    if start <= pos <= end:
                        key = (chrom, start, end)
                        window_meth[key] += mc
                        window_total[key] += total

        meth_series = pd.Series(window_meth, name=f"{sample_id}_meth")
        total_series = pd.Series(window_total, name=f"{sample_id}_depth")
        return sample_id, meth_series, total_series

    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return sample_id, pd.Series(dtype=int), pd.Series(dtype=int)

# Main parallel execution
sample_files = [os.path.join(allc_dir, f) for f in os.listdir(allc_dir) if f.endswith(".tsv.gz")]
sample_results = []

print("Processing ALLC files in parallel...")
with ProcessPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(process_allc_file, f) for f in sample_files]
    for future in tqdm(as_completed(futures), total=len(futures)):
        sample_id, meth_series, total_series = future.result()
        sample_results.append((sample_id, meth_series, total_series))

# Merge all results with reference
ref_reset['window_id'] = list(zip(ref_reset['chr'], ref_reset['start'], ref_reset['end']))
ref_reset = ref_reset.set_index('window_id')

for sample_id, meth_series, total_series in sample_results:
    ref_reset = ref_reset.join(meth_series, how='left')
    ref_reset = ref_reset.join(total_series, how='left')

# Fill missing values with zeros
ref_reset = ref_reset.fillna(0)

# Write output
print(f"Writing output to {out_path}...")
ref_reset.reset_index(drop=True).to_csv(out_path, sep='\t', index=False)
print("Done.")
