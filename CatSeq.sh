#!/usr/bin/env bash

if [[ $1 == *bam ]]; then
 ~/software/bedtools-2.17.0/bin/bamToFastq  -i $1 -fq /dev/stdout
elif [[ $1 == *SRA ]]; then
	/net/eichler/vol5/home/mchaisso/software/sratoolkit.2.3.3-4-ubuntu64/bin/fastq-dump $1 --stdout
elif [[ $1 == *fastq.gz ]]; then
  zcat $1
fi
