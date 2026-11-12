#!/bin/bash

# Check if input BED file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input.bed> [upstream] [downstream]"
    exit 1
fi

input_bed="$1"
output_bed="${input_bed%.bed}_promoters.bed"
upstream=${2:-1000}  # Default 1000bp upstream
 downstream=${3:-100} # Default 100bp downstream

# Generate promoter regions
awk -v up="$upstream" -v down="$downstream" 'BEGIN {OFS="\t"} {
    if ($6 == "+") {
        print $1, $2-up, $2+down, $4, $5, $6;
    } else {
        print $1, $3-down, $3+up, $4, $5, $6;
    }
}' "$input_bed" > "$output_bed"

echo "Promoter regions written to $output_bed"
