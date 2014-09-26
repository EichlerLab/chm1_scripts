#!/usr/bin/env bash
bed=$1

sam=$fasta.sam
rmsk=RepeatMasked/$fasta.out
blocks=$fasta.block
localbed=$(basename $bed)
base=${localbed%.*}
fasta=$base.fasta
echo $PBS/PrintInsertionContext.py $bed /var/tmp/mchaisso/ucsc.hg19.fasta $bed /var/tmp/mchaisso/ucsc.hg19.fasta $fasta

~/software/blasr_2/cpp/alignment/bin/blasr $fasta /var/tmp/mchaisso/ucsc.hg19.fasta -bwt /var/tmp/mchaisso/ucsc.hg19.fasta.bwt -ctab /var/tmp/mchaisso/ucsc.hg19.fasta.ctab -m 1 -out $fasta.m1 



#$PBS/SamToBlocks.py $sam $blocks


