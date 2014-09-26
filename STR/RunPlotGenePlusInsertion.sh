#!/usr/bin/env sh

mkdir -p $1
grep "$1" ~/projects/RefSeq/hg19.gene_list.bed > $1/exons.bed
grep "$1\$" /net/eichler/vol5/home/mchaisso/projects/RefSeq/hg19.gene_interval.bed > $1/coords.bed
bedtools intersect -a /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/STR.bed -b $1/coords.bed | tr -d ":FULL" > $1/insertion.str.bed
bedtools intersect -b $1/coords.bed -a /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/STR/hg19.str.full.bed > $1/hg19.str.bed


Rscript /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/STR/PlotGenePlusInsertions.R $1/exons.bed $1/hg19.str.bed $1/insertion.str.bed $1 $1.pdf 1
