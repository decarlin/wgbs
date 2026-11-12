#!/bin/bash

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/pd_allc_list.tsv \
	 --output_path /data/combined_wgbs/pd_celfie.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions celfie /home/user/celfie/celfie_regions.bed \
	 --quantifiers celfie count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/ad_allc_list.tsv \
	 --output_path /data/combined_wgbs/ad_celfie.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions celfie /home/user/celfie/celfie_regions.bed \
	 --quantifiers celfie count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/dlb_allc_list.tsv \
	 --output_path /data/combined_wgbs/dlb_celfie.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions celfie /home/user/celfie/celfie_regions.bed \
	 --quantifiers celfie count CGN
