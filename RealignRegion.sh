#!/usr/bin/env bash
if [ $# -lt 2 ]; then
  echo "usage: RealignRegion.sh alignemnts.bam region target"
	echo " region is in chr:start-end format"
	echo " target is the sequence that was assembled, and should be near the center of the region."
  exit
fi

module load python/latest
module load numpy/latest
source ~mchaisso/pacbio.env/bin/activate
d=$(dirname "$1")
samtools view $1 $2 > subset.sam
~mchaisso/software/blasr_2/cpp/pbihdfutils/bin/samtobas subset.sam subset.bas.h5
#~mchaisso/projects/PacBioSequencing/scripts/GetReadsInRegion.sh $1 $2 | ~mchaisso/projects/PacBioSequencing/scripts/FindBasFiles.py --dir $d --findex > reads.findex
~mchaisso/software/stable_blasr/blasr subset.bas.h5 $3 -sam -out subset.sam -bestn 1 -clipping subread
rm -f reads.fasta.cmp.h5
~mchaisso/software/stable_blasr/samtoh5 subset.sam $3 subset.cmp.h5 
cmph5tools.py sort subset.cmp.h5 --deep
~mchaisso/software/stable_blasr/loadPulses subset.bas.h5 subset.cmp.h5 -metrics MergeQV,DeletionQV,InsertionQV,SubstitutionQV,SubstitutionTag,DeletionTag -byread
quiver -o reads.consensus.fasta --referenceFilename $3 subset.cmp.h5 
