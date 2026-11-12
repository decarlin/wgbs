cd /home/user/UXM_deconv

python export_uxm_pctU.py \
       /data/ADRC_fastq/ADRC_hisat/pat_files/*pat.gz \
       -a /home/user/projects/wgbs/data/reference/EckerFIXED.UM100.l4.hg38.tsv \
       -o /data/ADRC_fastq/ADRC_hisat/uxm_ecker_pctU_output \
       --threads 1 \
       --verbose
