#!/usr/bin/env bash

if [[ $# < 2 ]]; then
   echo "usage: all_alignments.bam region_file [slop] "
   exit
fi

if [ $# > 3 ]; then
  slop="--slop $3"
else
  slop=""
fi

rm -f commands.txt
i=1
for line in `~mchaisso/projects/PacBioSequencing/scripts/GetRegionFromBed.py $2 $slop`; do
		if [ -e rgn_$i/region/9-terminator/region.ctg.fasta ]; then 
				
				echo "cd rgn_$i; ~/projects/PacBioSequencing/scripts/RealignRegion.sh $1 $line region/9-terminator/region.ctg.fasta; cd .. " >> commands.txt
		fi
		i=$(($i+1))
done

mkdir -p jobs
~mchaisso/scripts/jobify commands.txt jobs quiv -jobsPerFile 100
~mchaisso/scripts/submit_jobs -mem 4 jobs/quiv*.sh

