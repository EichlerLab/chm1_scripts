#!/bin/bash

/bin/rm -f commands.txt
for file in `cat $1`; do
 localfile=$(basename $file)
 base=${localfile%.*}
 echo "zcat $file | $PBS/queryFasta $2 /dev/stdin $base.counts"  >> commands.txt
done
mkdir -p jobs
/bin/rm -f jobs/*.sh
jobify commands.txt jobs gtp -jobsPerFile 4 
submit_jobs -slots 4 jobs/*.sh
