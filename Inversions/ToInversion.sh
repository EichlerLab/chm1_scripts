#!/usr/bin/env bash

cat $1/region/region.ctg.consensus.fasta | rc.pl > $1/rc.fasta
rgn=`cat $1/region.txt`
samtools faidx /var/tmp/mchaisso/ucsc.hg19.fasta `~/scripts/slop.py $rgn 5000` > $1/ref.fasta 
blasr $1/rc.fasta $1/ref.fasta -forwardOnly  -minPctIdentity 90 > $1/rev_alignment.m1
~/scripts/m1utils.py $1/rev_alignment.m1 --bed --out $1/rev_alignment.bed
blasr $1/region/region.ctg.consensus.fasta $1/ref.fasta -bestn 1 -sam -out $1/align.sam
samtools faidx $1/ref.fasta
$PBS/PrintGaps.py $1/ref.fasta $1/align.sam --outFile $1/align.bed
chr=`cat $1/region.txt | cut -f 1 -d":"`
start=`cat $1/region.txt | cut -f 2 -d":" | cut -f 1 -d'-'`
invStart=`cut -f 2 $1/rev_alignment.bed | head -1`
invEnd=`cut -f 3 $1/rev_alignment.bed | head -1`
echo $invEnd $invStart
echo $chr":"$((start+invStart))"-"$((start + invEnd))" "$((invEnd-invstart))
#echo "" | awk '{print $$chr":"$$start + $$invStart"-"$$start + $$invEnd" "$$invEnd - $$invStart }'
cat $1/align.bed | cut -f 1-5

