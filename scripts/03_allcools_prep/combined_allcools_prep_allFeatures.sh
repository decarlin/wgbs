#!/bin/bash

allcools generate-dataset  \
	 --allc_table /home/user/projects/wgbs/data/sample_lists/combined_file_list.tsv \
	 --output_path /data/combined_wgbs/outpc.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions chrom1k 1000 \
	 --regions eckerDMR /data/ecker_DMRs/dmr.MajorType.rho08.dist250.filter1.repflt_nochr.bed \
	 --regions genes /data/genome/gencode.v47.basic.annotation_coding_reordered_unique.bed \
	 --quantifiers genes count CGN \
	 --quantifiers eckerDMR count CGN \
	 --quantifiers chrom1k count CGN
