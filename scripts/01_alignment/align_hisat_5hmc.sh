#!/bin/bash
# cutadapt 
#ls /data/5hmc/*_L001_R1_001.fastq.gz | xargs -n1 basename | cut -d _ -f 1,2,3 > ~/5hmc_names.txt

export TMPDIR=/tmp
export tmp_dir=/tmp
export n_cores=16
export path_to_genome_folder=/data/genome/hisat_genome/hs2_genome
export picard=/home/user/picard.jar

while IFS= read -r line; do

    echo $line
    date

    echo "trimming..."
    
    trim_galore -j 8 --paired /data/5hmc/${line}_L001_R1_001.fastq.gz \
		/data/5hmc/${line}_L001_R2_001.fastq.gz \
		--three_prime_clip_R1 10 --clip_R2 10 -o /data/5hmc/tmp/

    #/data/full_align/2544-17_S9_L001_R1_001_val_1.fq.gz /data/full_align/2544-17_S9_L001_R2_001_val_2.fq.gz
   
    #just align to phage
    #methylpy build-reference \
        #       --input-files /data/lambda.fa.fasta  \
        #       --output-prefix /data/control_reference/lambda \
        #       --aligner bowtie2
    
    echo $line >> ~/5hmc_align_hisat_err.txt

    export bam_file=/data/5hmc/hisat/${line}.hisat3n_sorted.bam
    export bam_file_d=/data/5hmc/hisat/${line}.hisat3n_sorted_deduped.bam
    export dedup_summary=/data/5hmc/hisat/${line}_deduped_matrix.txt

    echo "mapping..."
    hisat-3n -x $path_to_genome_folder \
                -q \
                -1 /data/5hmc/tmp/${line}_L001_R1_001_val_1.fq.gz \
                -2 /data/5hmc/tmp/${line}_L001_R2_001_val_2.fq.gz \
                --base-change C,T \
                --no-repeat-index \
                --no-spliced-alignment \
                --directional-mapping \
                --no-temp-splicesite \
                -t -p $n_cores \
                --new-summary \
                --summary-file /data/5hmc/hisat/${line}_summary.txt | samtools sort -@ $n_cores -O BAM -o ${bam_file}

    echo "deduping..."
    java -Xmx32g -jar $picard MarkDuplicates \
	 I=${bam_file} \
	 O=${bam_file_d} \
	 M=${dedup_summary} \
	 VALIDATION_STRINGENCY=SILENT\
	 TMP_DIR=$temp_dir

    echo "get allc..."
    methylpy call-methylation-state \
                --input-file ${bam_file_d} \
                --path-to-output /data/5hmc/hisat/ \
                --paired-end True \
                --min-mapq 20 \
                --sample ${line} \
		--num-procs $n_cores \
               --ref-fasta /data/genome/hg38_lambda.fa

    
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R1_val_1.fq.gz
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R2_val_2.fq.gz

     echo $line |  mailx -s "5hmc hisat sample finished on cello-neuro ${line}" user@localhost
done < ~/5hmc_names.txt

echo "Done"
date
