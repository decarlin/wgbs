import os
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from tqdm import tqdm

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
ref['window_id'] = ref['chr'].astype(str) + ":" + ref['start'].astype(str) + "-" + ref['end'].astype(str)
ref = ref.set_index('window_id')

# Build interval lookup by chromosome
window_intervals = defaultdict(list)
for idx, row in ref.reset_index().iterrows():
    chrom = row['chr']
    start = row['start']
    end = row['end']
    window_id = row['window_id']
    window_intervals[chrom].append((start, end, window_id))

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

                for start, end, window_id in window_intervals[chrom]:
                    if start <= pos <= end:
                        window_meth[window_id] += mc
                        window_total[window_id] += total

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
with ProcessPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(process_allc_file, f) for f in sample_files]
    for future in tqdm(as_completed(futures), total=len(futures)):
        sample_id, meth_series, total_series = future.result()
        sample_results.append((sample_id, meth_series, total_series))

# Merge all results with reference
data = ref.copy()
for sample_id, meth_series, total_series in sample_results:
    data = data.join(meth_series, how='left')
    data = data.join(total_series, how='left')

# Fill missing values with zeros
data = data.fillna(0)

# Write output
print(f"Writing output to {out_path}...")
data.reset_index(drop=True).to_csv(out_path, sep='\t', index=False)
print("Done.")
