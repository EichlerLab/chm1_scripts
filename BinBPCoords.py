#!/usr/bin/env python

import argparse
import sys
import numpy as np
import Tools

ap = argparse.ArgumentParser(description="Bin bp coordinates.")
ap.add_argument("coords", help="Bed file of coordinates.")
ap.add_argument("fai", help="index of genome.")
ap.add_argument("col", help="which column to cluster [1 or 2]", type=int)
ap.add_argument("bin", help="Bin size", type=int)
ap.add_argument("--min", help="Min bin count to print (2)", type=int, default=2)
ap.add_argument("--ucsc", help="Print coordinates in UCSC browser compatible format.", action='store_true', default=False),
ap.add_argument("--expand", help="Widen the coordiantes by this.", type=int, default=0)
args = ap.parse_args()


fai = Tools.ReadFAIFile(args.fai)

def GetNumBins(length, binSize):
    if (length % binSize == 0):
        return length/binSize
    else:
        return length/binSize + 1

bins = dict( (chrom, np.zeros(GetNumBins(fai[chrom][0], args.bin))) for chrom in fai.keys() )
coordsFile = open(args.coords, 'r')
lineNumber = 0
for line in coordsFile:
    vals = line.split()
    chrom = vals[0]
    pos = int(vals[args.col])
    bins[chrom][pos/args.bin] += 1
    if (lineNumber % 100000 == 0):
        sys.stderr.write(str(lineNumber) + "\n")
    lineNumber += 1
    

    
for chrom in bins.keys():
    for i in range(len(bins[chrom])):
        if (bins[chrom][i] > args.min):
            start = max(0,args.bin*i-args.expand)
            end = min(fai[chrom][0], args.bin*(i+1)+args.expand)
            if (args.ucsc == False):
                output = "{}\t{}\t{}\t{}".format(chrom, start, end, bins[chrom][i])
            else:
                output = "{}:{}-{}\t{}".format(chrom, start, end, bins[chrom][i])
            print output
