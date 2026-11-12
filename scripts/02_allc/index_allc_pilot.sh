#ls -d /data/wgbs_pilot_2021/allc_*.tsv.gz >~/pilot_allc_list.txt

while IFS= read -r line; do

    echo $line
    
    gunzip -c "$line" > "${line%.gz}"
    bgzip --force ${line%.gz}
    allcools tabix-allc ${line}

done < ~/pilot_allc_list.txt
