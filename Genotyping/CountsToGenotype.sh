#!/usr/bin/env bash
for file in `ls *.counts`; do
base=${file%.*}
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Genotyping/CountsToGenotype.py $file $base.gt
echo $file
done
