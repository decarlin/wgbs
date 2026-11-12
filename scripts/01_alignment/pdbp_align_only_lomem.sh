#!/bin/sh

# cutadapt 
#ls /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/wgbs_072424/PDBP/*_R1_001.fastq.gz | xargs -n 1 basename | cut -d _ -f 1,2 | head -n 4> ~/wgbs_names_4.txt

#ls /data/pdbp_01302024/*_R1.fastq.gz |xargs -n 1 basename | cut -d _ -f 1 > ~/pdbp_names.txt

#ls /data/pdbp_01302024/trimmed/*_R1_val_1.fq.gz |xargs -n 1 basename | cut -d _ -f 1 >~/trimmed.txt

export TMPDIR=/data/tmp/

while IFS= read -r line; do

    echo $line
    date

    echo $line >> ~/pdbp_align_lomem_err.txt
    
    methylpy paired-end-pipeline \
         --trim-reads FALSE \
        --read1-files /data/pdbp_01302024/trimmed/${line}_R1_val_1.fq.gz \
        --read2-files /data/pdbp_01302024/trimmed/${line}_R2_val_2.fq.gz \
        --sample ${line} \
        --forward-ref /data/genome/combined_reference/hg38_lambda_f \
        --reverse-ref /data/genome/combined_reference/hg38_lambda_r \
        --ref-fasta /data/genome/hg38_lambda.fa \
        --num-procs 8 \
        --remove-clonal True \
        --aligner bowtie2 \
        --path-to-output /data/pdbp_01302024/aligned/ \
        --path-to-picard /home/user \
        --unmethylated-control chrL 2>> ~/pdbp_align_lomem_err.txt >>  ~/pdbp_align_lomem_stdout.txt

done < ~/pdbp_names_toLomem.txt

echo "Done"
date
