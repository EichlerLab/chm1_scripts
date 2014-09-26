#!/usr/bin/env sh
# first print gaps
$PBS/PrintGaps.py /var/tmp/mchaisso/ucsc.hg19.fasta $1 --tsd 10 --condense 20 --outFile $1.tmp
# second merge them
$PBS/rmdup.py $1.tmp $1.bed

rm $1.tmp
$PBS/GapBedToFasta.py $1.bed $1.insertion.fasta $1.deletion.fasta

mkdir -p $1.insertion
mkdir -p $1.deletion

RepeatMasker -dir $1.insertion -pa 8 $1.insertion.fasta
RepeatMasker -dir $1.deletion -pa 8 $1.deletion.fasta

grep "insertion" $1.bed > $1.insertion.bed
grep "deletion" $1.bed > $1.deletion.bed

$PBS/AnnotateGapBed.py $1.insertion.bed $1.insertion.annotated.bed $1.insertion/$1.insertion.fasta.out
$PBS/AnnotateGapBed.py $1.deletion.bed $1.deletion.annotated.bed $1.deletion/$1.deletion.fasta.out


