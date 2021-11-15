#!/bin/bash

# to do: clean python script
clusters=$1
fastq=$2
outdir=$3

mkdir -p $outdir
CLUSTERS_CNT=$(awk '($5 ~ /[0-9]/) {print $5}' $clusters | sort -nr | uniq | head -n1)
for ((i = 0 ; i <= $CLUSTERS_CNT ; i++));
do
    cluster_id=$i
    awk -v cluster="$cluster_id" '($5 == cluster) {print $1}' $clusters > $outdir/c$i.txt
    seqkit grep $fastq -f $outdir/c$i.txt -o $outdir/c$i.fastq
done