
ls /data/pdbp_01302024/aligned_hisat/bam_files/*.bam >~/pdbp_bams.txt

while IFS= read -r file; do

wgbstools bam2pat $file \
	  --out_dir /data/pdbp_01302024/aligned_hisat/pat_files

done < ~/pdbp_bams.txt
