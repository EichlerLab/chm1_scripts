#!/usr/bin/env python

import argparse
import Tools
import sys
import numpy as np

ap = argparse.ArgumentParser(description="Use a superdup file to compute maximal identity of a given position.")
ap.add_argument("bed", help="Input bed/duperdup file.")
ap.add_argument("genome", help="Genoem file with .fai")
ap.add_argument("--bin", help="Bin size", type=int, default=100)
ap.add_argument("--out", help="Output file.", default="/dev/null")
args = ap.parse_args()

bed = open(args.bed)

fai = Tools.ReadFAIFile(args.genome + ".fai")

outFile = open(args.out, 'w')

bins = {}
for k in fai.keys():
    l = fai[k][0]
    bins[k] = np.zeros(l/ args.bin + 1)

lineNumber = 0
for line in bed:
    vals = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end   = int(vals[2])
    ident = float(vals[25])
    if (chrom in bins):
        b = bins[chrom]
        startBin = start / args.bin
        endBin   = end   / args.bin
        
        for i in range(startBin, endBin):
            b[i] = max(b[i], ident)
    lineNumber += 1
    if (lineNumber % 100000 == 0):
        sys.stderr.write("processed {} \n".format(lineNumber))

for k in bins.keys():
    b = bins[k]
    print "writing " + str(len(b)) + " entries"
    i = 0
    rgnStart = i
    rgnIdentity = b[i]
    
    while (i < len(b)):
        if (i < len(b)-1 and b[i+1] == rgnIdentity):
            i+= 1
        else:
            if (b[rgnStart] > 0 ):
                outFile.write("{}\t{}\t{}\t{}\n".format(k, rgnStart*args.bin, (i+1)*args.bin, rgnIdentity))
            i+=1
            if (i < len(b)):
                rgnStart = i
                rgnIdentity = b[i]

outFile.close()
