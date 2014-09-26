#!/usr/bin/env sh
$PBS/GapBedToFasta.py $1 $1.fasta
~/software/bin/trf $1.fasta  2 7 7 80 10 20 500 -m -ngs -h  > $1.trf
b=${1%.*}
$PBS/AnnotateWithTRF.py $1 $1.trf $b.trf.bed
