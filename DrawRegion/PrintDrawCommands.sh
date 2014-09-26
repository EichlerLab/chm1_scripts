#!/usr/bin/env bash


echo $1
echo $2
echo $3
index=0
mkdir -p bedfiles


while read line; do
		query=`echo $line | tr -d '\n' `
		grep $query $1  | awk '{if ($9 == "C") { print $6"\t"$7"\t"0"\t"$10}  else { print  $6"\t"$7"\t"1"\t"$10}  }' > bedfiles/repeats.$index.bed
		start=`echo $query | cut -f 2 -d'/'`
		end=`echo $query | cut -f 3 -d'/'`
		chrom=`echo $query | cut -f 1 -d'/'`
		insert=$(($end - $start))
		echo $line
		# find the alignment in the blocks file
		grep $query $2 | awk '{print $3"\t"$4"\t"$5}' > bedfiles/pt.$index.bed
		# print repeat annotations
		egrep "\)n|Low|rich" bedfiles/repeats.$index.bed > bedfiles/str.$index.bed
		egrep -v "\)n|Low|rich" bedfiles/repeats.$index.bed > bedfiles/mei.$index.bed
		# draw the region
		mkdir -p plots
		echo Rscript $PBS/DrawRegion/DrawRegion.R --mei bedfiles/mei.$index.bed --str bedfiles/str.$index.bed --insert $insert --pos $start --window $3  --out $index.pdf --chrom $chrom --pt bedfiles/pt.$index.bed
		Rscript $PBS/DrawRegion/DrawRegion.R --mei bedfiles/mei.$index.bed --str bedfiles/str.$index.bed --insert $insert --pos $start --window $3  --out plots/$index.pdf --chrom $chrom --pt bedfiles/pt.$index.bed
		index=$((index + 1))

done 



