
ls -d /data/cfDNA_Jan2025/allc_*.tsv.gz >~/allc_list.txt

while IFS= read -r line; do

    echo $line
    
    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/allc_list.txt
