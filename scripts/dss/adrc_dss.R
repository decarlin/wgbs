library(arrow)
library(DSS)

meta<-read.csv('/data/combined_wgbs/combined_metadata.csv', header=TRUE, row.names=1)

sample_names <- rownames(meta)[meta$study=='ADRC']

feather_dir <- "/data/combined_wgbs/adrc_feather"

expected_files <- paste0(
    feather_dir,
    "/",
    sample_names,
    ".feather"
)

exists_vec <- file.exists(expected_files)

missing_samples <- sample_names[!exists_vec]

sample_names<-sample_names[!sample_names %in% missing_samples]

dat_list <- lapply(sample_names, function(s) {

    message("Loading ", s)

    df <- read_feather(
        paste0("/data/combined_wgbs/adrc_feather/", s, ".feather")
    )

    as.data.frame(df)
})

BSobj <- makeBSseqData(dat_list, sample_names)

design <- data.frame(
    dx = factor(meta$dx),
    sex = factor(meta$sex),
    study = factor(meta$study)
)