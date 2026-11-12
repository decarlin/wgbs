

while IFS= read -r samp; do

    export line=/data/ADRC_fastq/ADRC_hisat/allc_${samp}.tsv.gz
    
    cp $line /data/ADRC_fastq/ADRC_hisat/new_index/

done < ~/P117_redo_names.txt
