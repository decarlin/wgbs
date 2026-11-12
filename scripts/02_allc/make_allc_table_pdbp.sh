#!/bin/bash

# Check if directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

dir="$1"
output_file="pdbp_file_list.tsv"

#echo -e "Filename_Prefix\tFull_Path" > "$output_file"

# Loop through matching files
for file in "$dir"/allc_*.tsv.gz; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        prefix=${filename#allc_}
        prefix=${prefix%.tsv.gz}
        echo -e "PDBP$prefix\t$file" >> "$output_file"
    fi
done

echo "File list written to $output_file"
