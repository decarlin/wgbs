#!/bin/bash

# cutadapt
#pull from cold storage
#ls /media/user/Expansion1/PDBP_cold_storage/*_R1.fastq.gz | xargs -n 1 basename | cut -d _ -f 1 > ~/pdbp_cs_names.txt
# ls /media/user/Expansion1/PDBP_cold_storage/*_R1.fastq.gz | xargs -n 1 basename | cut -d _ -f 1 | head -n 1 > ~/pdbp_cs_names_TEST.txt

export TMPDIR=/tmp
export tmp_dir=/tmp
export n_cores=16
export path_to_genome_folder=/data/genome/hisat_genome/hs2_genome
export picard=~/picard.jar

while IFS= read -r line; do

    echo $line
    date

    echo "trimming..."
    
    trim_galore -j 8 --paired /media/user/Expansion1/PDBP_cold_storage/${line}_R1.fastq.gz \
    	/media/user/Expansion1/PDBP_cold_storage/${line}_R2.fastq.gz \
       --three_prime_clip_R1 10 --clip_R2 10 -o /data/pdbp_01302024/trimmed/

    #just align to phage
    #methylpy build-reference \
        #       --input-files /data/lambda.fa.fasta  \
        #       --output-prefix /data/control_reference/lambda \
        #       --aligner bowtie2
    
    echo $line >> ~/pdbp_align_hisat_err.txt

    export bam_file=/data/pdbp_01302024/aligned_hisat/bam_files/${line}.hisat3n_sorted.bam
    export bam_file_d=/data/pdbp_01302024/aligned_hisat/bam_files/${line}.hisat3n_sorted_deduped.bam
    export dedup_summary=/data/pdbp_01302024/aligned_hisat/bam_files/${line}_deduped_matrix.txt

    echo "mapping..."
    hisat-3n -x $path_to_genome_folder \
                -q \
                -1 /data/pdbp_01302024/trimmed/${line}_R1_val_1.fq.gz \
                -2 /data/pdbp_01302024/trimmed/${line}_R2_val_2.fq.gz \
                --base-change C,T \
                --no-repeat-index \
                --no-spliced-alignment \
                --directional-mapping \
                --no-temp-splicesite \
                -t -p $n_cores \
                --new-summary \
                --summary-file /data/pdbp_01302024/aligned_hisat/${line}_summary.txt | samtools sort -@ $n_cores -O BAM -o ${bam_file}

    echo "deduping..."
    java -Xmx32g -jar $picard MarkDuplicates \
	 I=${bam_file} \
	 O=${bam_file_d} \
	 M=${dedup_summary} \
	 TMP_DIR=$temp_dir \
	 ASSUME_SORTED=true \
	 VALIDATION_STRINGENCY=SILENT

    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R1_val_1.fq.gz
    #rm /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/tmp/${line}_R2_val_2.fq.gz

     [ -n "$NOTIFY_EMAIL" ] && echo $line | mailx -s "hisat sample finished on cello-neuro" "$NOTIFY_EMAIL"
done < ~/pdbp_cs_names.txt

echo "Done"
date
