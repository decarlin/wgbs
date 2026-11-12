# WGBS cfDNA methylation analysis

Scripts for processing whole-genome bisulfite sequencing (WGBS) cell-free DNA from the ADRC, PDBP, and pilot cohorts. They cover alignment, ALLC generation, ALLCools feature matrices, PCA and classifier analyses (AD, PD, DLB vs. controls), and cell-type deconvolution with CelFiE and UXM.

> Scripts use absolute paths on the analysis server (`/data/...`, `~/...`, `/home/user/projects/wgbs/data/...`). They are a record of the analysis, not a portable pipeline.
>
> Some alignment scripts email a notification after each sample when `NOTIFY_EMAIL` is set (e.g. `export NOTIFY_EMAIL=you@example.org`); otherwise they skip it.

## Layout

| Directory | Contents |
|---|---|
| `scripts/01_alignment/` | Trimming (trim_galore) and hisat-3n alignment; BAM generation and cleanup |
| `scripts/02_allc/` | Building and indexing ALLC tables |
| `scripts/03_allcools_prep/` | `allcools generate-dataset` runs (1k/100k bins, CHN, gene regions, Ecker DMRs) and file-list/metadata prep |
| `scripts/04_allcools_features/` | Loading MCDS into AnnData, QC filtering, feature matrices per cohort |
| `scripts/05_analysis/` | PCA, covariate models, classifiers/ROC, longitudinal analyses, heatmaps, boxplots |
| `scripts/celfie/` | CelFiE input prep and run scripts |
| `scripts/uxm/` | UXM / wgbstools (`bam2pat`) deconvolution |
| `scripts/dss/` | ALLC → DSS conversion and DSS differential methylation |
| `scripts/utils/` | Small helpers (promoters BED, chrom sizes, unique BED) |

Not tracked in git (see `.gitignore`):

- `data/sample_lists/`: ALLC file lists per cohort (sample IDs → server paths)
- `data/reference/`: Ecker UXM marker atlas and per-region/per-cell-type DMR tables
- `data/results/`: assorted outputs
- `archive/`: superseded script versions, editor backups and autosaves, and junk files. `archive/MOVES.tsv` records where every file came from.

## Naming

Dated scripts use an `_MMDDYYYY` suffix. Only the latest version of each script is kept in `scripts/`; older versions are in `archive/superseded/` locally and will be tracked by git from now on.

## Dependencies

- Python: pandas, numpy, scanpy, ALLCools, scikit-learn, statsmodels, pybedtools, matplotlib, seaborn, xarray, diffxpy, tqdm
- R: DSS, arrow, stringr
- CLI: trim_galore, hisat-3n, allcools, wgbstools, CelFiE (`~/celfie/scripts/celfie.py`), Java
