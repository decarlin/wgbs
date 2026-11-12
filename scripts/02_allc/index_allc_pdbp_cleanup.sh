#ls -d /data/pdbp_01302024/aligned_hisat/allc_*.tsv.gz >~/pdbp_allc_list.txt

while IFS= read -r line; do

    echo $line
    
    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/pdbp_allc_list.txt
