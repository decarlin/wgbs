

ls /data/ADRC_fastq/ADRC_hisat/bam_files/*.bam >~/adrc_bams.txt

while IFS= read -r file; do

wgbstools bam2pat $file \
	  --out_dir /data/ADRC_fastq/ADRC_hisat/pat_files

done < ~/adrc_bams.txt
