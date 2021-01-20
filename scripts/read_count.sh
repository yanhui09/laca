#!/bin/bash
# shell script to summarize classified reads
        
awk 'BEGIN {FS="\t"}; {print $2}' "$1/adjusted-$2_tax_assignments.txt" > "$1/taxa.tmp"
num=$(wc -l "$1/taxa.tmp" | awk '{print $1}')-1
for ((i=0; i<=$num; i++)); do echo 1; done > "$1/col2"
for ((i=0; i<=$num; i++)); do echo ONTU_$i; done > "$1/col1"
        
paste "$1/col1" "$1/col2" "$1/taxa.tmp" > "$1/step1.tmp"
        
TAB=$'\t'
HASH=$'#'
echo "$HASH OTU_ID $TAB $2 $TAB taxonomy" | sed "s/ //g" > "$1/header.tmp"
cat "$1/header.tmp" "$1/step1.tmp" | sed "s/^$/Unassigned/g"> "$1/$2.tmp"

#clean
rm -f "$1/col1" "$1/col2" "$1/taxa.tmp" "$1/header.tmp" "$1/step1.tmp"        