#!/usr/bin/env bash

echo "" > commands.txt

for file in `cat /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/Genotyping/pcr_free.txt`; do
		a=$(basename $file)
		b=${a%.*}
		echo "bamToFastq -i $file -fq /dev/stdout | $PBS/queryFasta $1 /dev/stdin $b.counts" >> commands.txt
done


echo "$PBS/queryFasta $1 /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/Validation/Venter/venter.fasta  Venter.counts" >> commands.txt

for file in `ls /net/eichler/vol19/projects/CHM1_project/nobackups/UWIllumina/s_*.fq.gz`; do
		a=$(basename $file)
		b=${a%.*}
		echo "zcat $file | $PBS/queryFasta $1 /dev/stdin $b.counts" >> commands.txt
done
