

while IFS= read -r samp; do

    export line=/data/ADRC_fastq/ADRC_hisat/new_index/allc_${samp}.tsv.gz
    
    echo $line

    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/P117_redo_names.txt
