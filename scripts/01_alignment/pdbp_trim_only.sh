
# cutadapt 
#ls /datasets/Public_Datasets/Desplats_cfDNA_Jun2024/wgbs_072424/PDBP/*_R1_001.fastq.gz | xargs -n 1 basename | cut -d _ -f 1,2 | head -n 4> ~/wgbs_names_4.txt

#ls /data/pdbp_01302024/*_R1.fastq.gz |xargs -n 1 basename | cut -d _ -f 1 > ~/pdbp_names.txt

while IFS= read -r line; do

    echo $line
    date
    
    trim_galore -j 4 --paired /data/pdbp_01302024/${line}_R1.fastq.gz \
		/data/pdbp_01302024/${line}_R2.fastq.gz \
		--three_prime_clip_R1 10 --clip_R2 10 -o /data/pdbp_01302024/trimmed

    #/data/full_align/2544-17_S9_L001_R1_001_val_1.fq.gz /data/full_align/2544-17_S9_L001_R2_001_val_2.fq.gz
   
    #just align to phage
    #methylpy build-reference \
        #       --input-files /data/lambda.fa.fasta  \
        #       --output-prefix /data/control_reference/lambda \
        #       --aligner bowtie2

done < ~/pdbp_names_toTrim.txt

echo "Done"
date
