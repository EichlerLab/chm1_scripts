#!/usr/bin/env bash
i=0
while [ $i -lt $3 ]; do
bedtools shuffle -i $1 -g /var/tmp/mchaisso/ucsc.hg19.fasta.fai  | bedtools sort | bedtools intersect -a stdin -b $2 -u  | wc -l
i=$((i+1))
done

