
ls -d /data/ADRC_fastq/ADRC_hisat/allc_*.tsv.gz >~/adrc_allc_list.txt

while IFS= read -r line; do

    echo $line
    
    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/adrc_allc_list.txt
