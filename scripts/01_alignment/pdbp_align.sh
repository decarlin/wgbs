
# cutadapt 
#ls /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/wgbs_072424/PDBP/*_R1_001.fastq.gz | xargs -n 1 basename | cut -d _ -f 1,2 | head -n 4> ~/wgbs_names_4.txt

ls /data/pdbp_01302024/*_R1.fastq.gz |xargs -n 1 basename | cut -d _ -f 1 > ~/pdbp_names.txt

while IFS= read -r line; do

    echo $line
    date
    
    trim_galore -j 8 --paired /data/pdbp_01302024/${line}_R1.fastq.gz \
		/data/pdbp_01302024/${line}_R2.fastq.gz \
		--three_prime_clip_R1 10 --clip_R2 10 -o /home/user/tmp/

    #/data/full_align/2544-17_S9_L001_R1_001_val_1.fq.gz /data/full_align/2544-17_S9_L001_R2_001_val_2.fq.gz
   
    #just align to phage
    #methylpy build-reference \
        #       --input-files /data/lambda.fa.fasta  \
        #       --output-prefix /data/control_reference/lambda \
        #       --aligner bowtie2

    methylpy paired-end-pipeline \
         --trim-reads FALSE \
        --read1-files /home/user/tmp/${line}_R1_val_1.fq.gz \
        --read2-files /home/user/tmp/${line}_R2_val_2.fq.gz \
        --sample ${line} \
        --forward-ref /data/genome/combined_reference/hg38_lambda_f \
        --reverse-ref /data/genome/combined_reference/hg38_lambda_r \
        --ref-fasta /data/genome/hg38_lambda.fa \
        --num-procs 16 \
        --remove-clonal True \
        --aligner bowtie2 \
        --path-to-output /data/pdbp_01302024/aligned/ \
        --path-to-picard /home/user \
        --unmethylated-control chrL 2> ~/adrc_align_err.txt >>  ~/adrc_align_stdout.txt

    rm /home/user/tmp/${line}_R1_val_1.fq.gz
    rm /home/user/tmp/${line}_R2_val_2.fq.gz

done < ~/pdbp_names.txt

echo "Done"
date
