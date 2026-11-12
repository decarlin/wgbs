#!/bin/bash
# cutadapt 
#ls /data/ADRC_fastq/ADRC_Pilot/*_R1.fastq.gz | xargs -n 1 basename | cut -d _ -f 1,2 > ~/Pilot_names.txt

export TMPDIR=/tmp
export tmp_dir=/tmp
export n_cores=28
export path_to_genome_folder=/data/genome/hisat_genome/hs2_genome
export picard=/home/user/picard.jar

while IFS= read -r line; do

    echo $line
    date

    echo "trimming..."
    
    trim_galore -j 8 --paired /data/ADRC_fastq/ADRC_Pilot/${line}_R1.fastq.gz \
		/data/ADRC_fastq/ADRC_Pilot/${line}_R2.fastq.gz \
		--three_prime_clip_R1 10 --clip_R2 10 -o /data/ADRC_fastq/tmp/

    #/data/full_align/2544-17_S9_L001_R1_001_val_1.fq.gz /data/full_align/2544-17_S9_L001_R2_001_val_2.fq.gz
   
    #just align to phage
    #methylpy build-reference \
        #       --input-files /data/lambda.fa.fasta  \
        #       --output-prefix /data/control_reference/lambda \
        #       --aligner bowtie2
    
    echo $line >> ~/adrc_align_hisat_err.txt

    export bam_file=/data/ADRC_fastq/ADRC_hisat/${line}.hisat3n_sorted.bam
    export bam_file_d=/data/ADRC_fastq/ADRC_hisat/${line}.hisat3n_sorted_deduped.bam
    export dedup_summary=/data/ADRC_fastq/ADRC_hisat/${line}_deduped_matrix.txt

    echo "mapping..."
    hisat-3n -x $path_to_genome_folder \
                -q \
                -1 /data/ADRC_fastq/tmp/${line}_R1_val_1.fq.gz \
                -2 /data/ADRC_fastq/tmp/${line}_R2_val_2.fq.gz \
                --base-change C,T \
                --no-repeat-index \
                --no-spliced-alignment \
                --directional-mapping \
                --no-temp-splicesite \
                -t -p $n_cores \
                --new-summary \
                --summary-file /data/ADRC_fastq/ADRC_hisat/${line}_summary.txt | samtools sort -@ $n_cores -O BAM -o ${bam_file}

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
                --path-to-output /data/ADRC_fastq/ADRC_hisat/ \
                --paired-end True \
                --min-mapq 20 \
                --sample ${line} \
		--num-procs $n_cores \
               --ref-fasta /data/genome/hg38_lambda.fa

    
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R1_val_1.fq.gz
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R2_val_2.fq.gz

     echo $line |  mailx -s "hisat sample finished on cello-neuro" user@localhost
done < ~/Pilot_names.txt

echo "Done"
date
