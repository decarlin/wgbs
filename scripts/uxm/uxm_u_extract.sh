cd /home/user/UXM_deconv

python export_uxm_pctU.py \
       /data/ADRC_fastq/ADRC_hisat/pat_files/*pat.gz \
       -a /home/user/UXM_deconv/supplemental/Atlas.U25.l4.hg38.blood_neuron_oligodend.tsv \
       -o /data/ADRC_fastq/ADRC_hisat/uxm_pctU_output \
       --threads 1 \
       --verbose
