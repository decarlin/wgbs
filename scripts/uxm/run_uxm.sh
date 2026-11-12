
uxm deconv /data/ADRC_fastq/ADRC_hisat/pat_files/*pat.gz -o /data/ADRC_fastq/ADRC_hisat/uxm_output.csv -a /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.tsv

uxm plot /data/ADRC_fastq/ADRC_hisat/uxm_output.csv -o /data/ADRC_fastq/ADRC_hisat/adrc_uxm_output.pdf
