#!/bin/bash

allcools generate-dataset  \
#--allc_table /home/user/projects/wgbs/data/sample_lists/file_list.tsv \
--output_path /data/cfDNA_Jan2025/allcools/outpc.mcds \
--chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
--obs_dim cell  \
--cpu 16 \
--chunk_size 50 \
--regions chrom1k 1000 \
--quantifiers chrom1k count CGN
