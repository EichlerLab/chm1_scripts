#!/usr/bin/env csh


samtools view $1 $2 | cut -f 1 | awk -F"/" '{print $1"/"$2}'  | sort | uniq
