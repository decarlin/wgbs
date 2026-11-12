#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input.fa>"
    exit 1
fi

input_file="$1"
output_file="${input_file%.fa}.sizes"

echo "Generating chromosome sizes..."
awk '/^>/ {if (seqlen){print chr,seqlen} chr=$0; sub(/^>/, "", chr); seqlen=0; next} {seqlen += length($0)} END {print chr,seqlen}' "$input_file" > "$output_file"

echo "Chromosome sizes written to $output_file"
