#!/usr/bin/env bash


file=$1
extension="${file##*.}"

if [ $extension = "gz" ]; then
		zcat $file
elif [ $extension = "bam" ]; then
		bamToFastq -i $file -fq /dev/stdout
elif [ $extension = "fastq" ]; then
		cat $file
fi
