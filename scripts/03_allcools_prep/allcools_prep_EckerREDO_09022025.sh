#!/bin/bash

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/pd_allc_list.tsv \
	 --output_path /data/combined_wgbs/pd_ecker.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 20 \
	 --chunk_size 20 \
	 --regions eckerDMR /data/ecker_DMRs/dmr.MajorType.rho08.dist250.filter1.repflt_nochr.bed.gz \
	 --quantifiers eckerDMR count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/ad_allc_list.tsv \
	 --output_path /data/combined_wgbs/ad_ecker.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 20 \
	 --chunk_size 20 \
	 --regions eckerDMR /data/ecker_DMRs/dmr.MajorType.rho08.dist250.filter1.repflt_nochr.bed.gz \
	 --quantifiers eckerDMR count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/dlb_allc_list.tsv \
	 --output_path /data/combined_wgbs/dlb_ecker.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 20 \
	 --chunk_size 20 \
	 --regions eckerDMR /data/ecker_DMRs/dmr.MajorType.rho08.dist250.filter1.repflt_nochr.bed.gz \
	 --quantifiers eckerDMR count CGN
