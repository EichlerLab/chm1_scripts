#!/usr/bin/env bash

for dir in `ls -d rgn_*`; do
 echo $dir
 rgn=`cat $dir/region.txt`;
 basename=`$HOME/scripts/rgn2name $rgn`
 $PBS/DotPlot.py --query $dir/region/region.ctg.consensus.fasta --target /var/tmp/mchaisso/ucsc.hg19.fasta --region $rgn --slop 2000 --savefig $basename.png --nolegend --matches dot:13 & 
done
