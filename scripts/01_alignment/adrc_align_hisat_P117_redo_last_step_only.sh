#!/bin/bash
# cutadapt 
#ls /data/ADRC_fastq/ADRC_P117/*_R1.fastq.gz | xargs -n 1 basename | cut -d _ -f 1 > ~/P117_names.txt

export TMPDIR=/tmp
export tmp_dir=/tmp
export n_cores=28
export path_to_genome_folder=/data/genome/hisat_genome/hs2_genome
export picard=/home/user/picard.jar

while IFS= read -r line; do

    echo $line
    date

    export bam_file_d=/data/ADRC_fastq/ADRC_hisat/${line}.hisat3n_sorted_deduped.bam
 
    echo "get allc..."
    methylpy call-methylation-state \
                --input-file ${bam_file_d} \
                --path-to-output /data/ADRC_fastq/ADRC_hisat/ \
                --paired-end True \
                --min-mapq 20 \
                --sample ${line} \
		--num-procs $n_cores \
               --ref-fasta /data/genome/hg38_lambda.fa

    
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R1_val_1.fq.gz
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R2_val_2.fq.gz

     echo $line |  mailx -s "hisat sample finished on cello-neuro ${line}" user@localhost
done < ~/P117_redo_names_AGAIN.txt

echo "Done"
date
