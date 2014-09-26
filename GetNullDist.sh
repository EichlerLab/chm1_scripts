#!/usr/bin/env python

bedtools shuffle -i $1 -g /var/tmp/mchaisso/ucsc.hg19.fasta.fai  | bedtools sort | bedtools intersect -a stdin -b $2 -u  | wc -l
