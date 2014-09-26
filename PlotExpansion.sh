#!/usr/bin/env bash
echo "$1\$"
grep "$1\$" /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/STR/refseq_genes_ucsc_01082014.bed | head -1 > $1.gene.bed
bedtools intersect -a ~/projects/RefSeq/hg19.gene_list.merged.bed -b $1.gene.bed | bedtools sort > $1.exons.bed
#grep $1 ~/projects/ExomeAssembly/RefSeq/ExonList.full.bed  | bedtools sort | bedtools merge > $1.exons.bed 
bedtools intersect -a /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/STR/hg19.str.full.bed -b $1.gene.bed > $1.str.bed
bedtools intersect -a  /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/STR/STR.bed -b $1.gene.bed | cut -f 1-10 > $1.str_expand.bed 
bedtools intersect -a  /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/TRF.bed -b $1.gene.bed |cut -f 1-10 >> $1.str_expand.bed

echo Rscript $PBS/PlotGenePlusInsertions.R  $1.exons.bed $1.str.bed $1.str_expand.bed $1.pdf  1
Rscript $PBS/PlotGenePlusInsertions.R $1.exons.bed $1.str.bed $1.str_expand.bed $1 $1.pdf 1
