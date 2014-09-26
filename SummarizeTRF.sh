#!/usr/bin/env bash

i=`cut -f1-5 $1 | bedtools intersect -a stdin -b /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/interstitial.bed | cut -f 5 | stats.py   | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "interstitial $i"
c=`cut -f1-5 $1 | bedtools intersect -a stdin -b /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/centromeric.hg19.bed | cut -f 5 | stats.py   | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "centromeric $c"
t=`cut -f1-5 $1 | bedtools intersect -a stdin -b /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/telomeric.bed | cut -f 5 | stats.py   | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "telomeric $t"

