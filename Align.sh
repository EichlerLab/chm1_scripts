#!/usr/bin/env bash

~/software/blasr_2/cpp/alignment/bin/blasr $1 /var/tmp/mchaisso/ucsc.hg19.fasta -sa /var/tmp/mchaisso/ucsc.hg19.fasta.sa -ctab /var/tmp/mchaisso/ucsc.hg19.fasta.ctab -sam  -bestn 1 -out $1.sam
samtools view -bS $1.sam | samtools sort - $1
samtools index $1.bam
