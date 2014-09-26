#!/usr/bin/env python

import pysam
import sys
import argparse
import numpy as np

ap = argparse.ArgumentParser(description="Print the average mapping quality for a list of exons.")
ap.add_argument("bam", help="Input alignments.")
ap.add_argument("bed", help="Input exon list in BED format.")
ap.add_argument("-o", "--out", help="Output file.", default=None)
ap.add_argument("-m", "--min", help="Minimum mapqv to count.", default=5, type=int)
args = ap.parse_args()

if (args.out == None):
    outFile = sys.stdout
else:
    outFile = open(args.out, 'w')

bamfile = pysam.Samfile( args.bam, "rb" )

bedfile = open(args.bed)
for line in bedfile:
    vals = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end   = int(vals[2])
    mapqv = np.zeros(end-start)
    count = np.zeros(end-start)
    for col in bamfile.pileup( chrom, start, end):
        
        if (col.pos < start or col.pos >= end):
            continue
        for read in col.pileups:
            if (read.alignment.mapq > args.min):
                mapqv[col.pos - start] = 1

    mean = np.sum(mapqv) / float(end - start)
    outFile.write("\t".join(vals) + "\t" + "{:2.2f}".format(mean) + "\n")


if (args.out is not None):
    outFile.close()
    
