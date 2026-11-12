#!/bin/bash

allcools generate-dataset  \
--allc_table file_list.tsv \
--output_path /data//allcools/outpc_genes.mcds \
--chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
--obs_dim cell  \
--cpu 16 \
--chunk_size 50 \
--regions genes /data/genome/gencode.v47.basic.annotation_coding_reordered_unique.bed \
--quantifiers genes count CGN
