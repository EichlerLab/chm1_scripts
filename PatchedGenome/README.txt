1. Create a seq map.

The following creates the seq map to go from str insertion coordinates to original hg19 coordinates. 

 bedtools intersect -a /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/STR.bed -b  /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/STR/hg19.str.full.bed -loj  | $PBS/rmdup.py /dev/stdin /dev/stdout --leftjustify | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$11"\t"$12"\t"$13}' > STR.modmap.bed

For MEI, it si only necessary to create insertions, 

cat AluY.simple.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$2"\t"$2+1'} > AluY.simple.modmap.bed
