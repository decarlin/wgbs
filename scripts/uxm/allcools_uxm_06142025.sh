#!/bin/bash

tail -n +2 /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.tsv | sed 's/^chr//' > /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.bed


allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/pd_allc_list.tsv \
	 --output_path /data/combined_wgbs/pd_uxm.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions uxm /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.bed \
	 --quantifiers uxm count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/ad_allc_list.tsv \
	 --output_path /data/combined_wgbs/ad_uxm.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions uxm /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.bed \
	 --quantifiers uxm count CGN

allcools generate-dataset  \
	 --allc_table /data/combined_wgbs/dlb_allc_list.tsv \
	 --output_path /data/combined_wgbs/dlb_uxm.mcds \
	 --chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
	 --obs_dim cell  \
	 --cpu 16 \
	 --chunk_size 50 \
	 --regions uxm /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.bed \
	 --quantifiers uxm count CGN
