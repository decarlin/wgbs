#!/bin/bash

allcools generate-dataset  \
--allc_table /home/user/projects/wgbs/data/sample_lists/pdbp_file_list.tsv \
--output_path /data/pdbp_01302024/allcools/pdbp_eckerDMR.mcds \
--chrom_size_path /data/genome/hg38_lambda_normal_noCHR.sizes \
--obs_dim cell  \
--cpu 3 \
--chunk_size 50 \
--regions eckerDMR /data/ecker_DMRs/dmr.MajorType.rho08.dist250.filter1.repflt_nochr.bed \
--quantifiers eckerDMR count CGN
