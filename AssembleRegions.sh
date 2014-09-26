#!/usr/bin/env bash

i=1
/bin/rm -f commands.txt
if [[ $# < 1 ]]; then
		echo "usage AssembleRegions.sh regions.bed [slop] [slop outputdir]"
		exit
fi

if [[ $# > 1 ]]; then
  slop="--slop $2"
else
  slop=""
fi
if [[ $# > 2 ]]; then
  output=$3
else
  output=output
fi

rm -f commands.txt
for line in `~mchaisso/projects/PacBioSequencing/scripts/GetRegionFromBed.py $1 $slop`; do
		mkdir -p $output/rgn_$i/region
		echo $line > $output/rgn_$i/region.txt
		i=$(($i+1))
done

