#!/usr/bin/env bash


/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Genotyping/QueriesToFasta.py $1 $1.fasta 
mrfast --search /net/eichler/vol2/eee_shared/assemblies/hg19/indexes/mrfast2.5/hg19 --seq $1.fasta -e 1  -o $1.fasta.sam 
grep -v "@" $1.fasta.sam | awk '{ if ($6 == "25M") print;}' | cut -f 10 | uniq -c | awk '{ print $2"\t"$1 }' > $1.fasta.degenerate
