#!/bin/bash

input="/data/genome/gencode.v47.basic.annotation_coding_reordered.bed"
output="/data/genome/gencode.v47.basic.annotation_coding_reordered_unique.bed"

awk 'BEGIN{OFS="\t"} {
  name = $4
  count[name]++
  if (count[name] == 1) {
    print $0
  } else {
    $4 = name "_" count[name]
    print $0
  }
}' "$input" > "$output"
