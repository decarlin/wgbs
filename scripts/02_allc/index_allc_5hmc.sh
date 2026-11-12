
ls -d /data/5hmc/hisat/allc_*.tsv.gz >~/5hmc_allc_list.txt

while IFS= read -r line; do

    echo $line
    
    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/5hmc_allc_list.txt
