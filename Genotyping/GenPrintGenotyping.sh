#!/bin/bash


for file in `cat $1`; do
 localfile=$(basename $file)
 base=${localfile%.*}
 echo "bamToFastq -i $file -fq /dev/stdout | $PBS/queryFasta $2 /dev/stdin $base.counts" 
done
