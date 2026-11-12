#!/bin/bash

allcools generate-dataset  \
--allc_table file_list.tsv \
--output_path /data/cfDNA_Jan2025/allcools/outpc.mcds \
--chrom_size_path /data/genome/hg38_lambda_normal.sizes \
--obs_dim cell  \
--cpu 16 \
--chunk_size 50 \
--regions chrom100k 100000 \
--regions chrom1k 1000 \
--regions genes /data/genome/gencode.v47.basic.annotation.bed \
--regions promoters /data/genome/gencode.v47.basic.annotation_coding_promoters.bed \
--quantifiers chrom100k count CGN,CHN \
--quantifiers chrom1k count CGN \
--quantifiers genes count CGN \
--quantifiers promoters count CGN
