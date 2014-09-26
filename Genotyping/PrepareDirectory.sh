#!/usr/bin/env sh


$PBS/Genotyping/CombineCountFiles.py --input s_*.counts --output CHM1.UW.counts
mkdir -p chm1_backup
mv s_*.counts chm1_backup
for file in `ls *.counts`; do
		echo $file
		base=${file%.*}
		$PBS/Genotyping/CountsToGenotype.py $file $base.gt --max 200 
done


