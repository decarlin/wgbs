import os
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from tqdm import tqdm
from pybedtools import BedTool

tissues = [
    "dendritic_cells", "endothelial_cells", "eosinophils", "erythroblasts",
    "macrophages", "monocytes", "neutrophils", "placenta", "T_cells",
    "adipose", "brain", "fibroblasts", "heart_left_ventricle", "hepatocytes",
    "lung", "mammary_gland", "megakaryocytes", "skeletal_muscle_myoblasts", "small_intestine"
]

ref_path = "/home/user/celfie/TIMs/sample_tims_summed.txt"

# Load reference
print("Loading reference data...")
ref_cols = ['chr', 'start', 'end'] + [f"{t}_meth" for t in tissues] + [f"{t}_depth" for t in tissues]
ref = pd.read_csv(ref_path, sep='\t', header=None)
ref.columns = ref_cols
ref = ref.astype({"start": int, "end": int})

ref['name'] = ref[['chr', 'start', 'end']].astype(str).agg('_'.join, axis=1)
ref['chr']=ref['chr'].str.replace("^chr", "", regex=True)

ref[['chr', 'start', 'end','name']].to_csv('/home/user/celfie/celfie_regions.bed', sep='\t', header=False, index=False)

