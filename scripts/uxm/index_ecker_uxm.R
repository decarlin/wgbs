uxm_markers<-read.csv('/home/user/projects/wgbs/data/reference/Ecker.UM100.l4.hg38.tsv', sep='\t',header=TRUE)

cpg <- read.table(
    gzfile("/home/user/wgbs_tools/references/hg38/CpG.bed.gz"),
    sep = "\t",
    col.names = c("chr", "pos", "cpg_index")
)

cpg_by_chr <- split(cpg, cpg$chr)

uxm_markers$startCpG <- NA_integer_
uxm_markers$endCpG <- NA_integer_

for (chr in names(cpg_by_chr)) {

    idx <- which(uxm_markers$chr == chr)

    if (length(idx) == 0)
        next

    chr_cpg <- cpg_by_chr[[chr]]

    # First CpG >= start
    start_idx <- findInterval(
        uxm_markers$start[idx],
        chr_cpg$pos
    ) + 1

    # Last CpG <= end
    end_idx <- findInterval(
        uxm_markers$end[idx],
        chr_cpg$pos
    )

    # Keep only valid indices
    valid_start <- start_idx <= nrow(chr_cpg)
    valid_end <- end_idx >= 1

    uxm_markers$startCpG[idx[valid_start]] <-
        chr_cpg$cpg_index[start_idx[valid_start]]

    uxm_markers$endCpG[idx[valid_end]] <-
        chr_cpg$cpg_index[end_idx[valid_end]]
}

write.table(uxm_markers,file='/home/user/projects/wgbs/data/reference/EckerFIXED.UM100.l4.hg38.tsv',row.names=FALSE,quote=FALSE, sep='\t')
