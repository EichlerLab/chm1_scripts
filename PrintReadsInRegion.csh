#!/usr/bin/env bash
echo $#
if [ $# -lt 6 ]; then
  echo "usage PrintReadsInRegion.csh alignments.bam chrom start end outdir path_to_reads"
	exit
fi
mkdir -p $5
region=$2:$3-$4
for read in `samtools view -q 20  $1 $2:$3-$4  | cut -f 1`; do
  echo $2 "	" $3 "	" $4 "	" $read >> $5/$2.$3.$4.bed
done


~/projects/PacBioSequencing/scripts/RunDotPlots.py  $5/$2.$3.$4.bed /var/tmp/mchaisso/ucsc.hg19.fasta $5 $6  --expand 1500 --j 4
