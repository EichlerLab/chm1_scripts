#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Print fraction of reference or read matched.")
ap.add_argument("m1", help="Input m1 file.")
ap.add_argument("--frac", "-f", help="Minimum fraction aligned (of max sink/read)", type=float, default=0.90)

args = ap.parse_args()

m1File = open(args.m1)

for line in m1File:
    vals = line.split()
    readLen = int(vals[11])
    readAlnStart = int(vals[9])
    readAlnEnd   = int(vals[10])
    readAlnLen = float(readAlnEnd - readAlnStart)
    readAlnFrac = float(readAlnEnd - readAlnStart)/readLen

    
    sinkLen = int(vals[8])
    sinkAlnStart = int(vals[6])
    sinkAlnEnd = int(vals[7])
    sinkAlnFrac = float(sinkAlnEnd - sinkAlnStart)/sinkLen

    maxFrac = max(readAlnFrac, sinkAlnFrac)
    
    if (maxFrac == readAlnFrac):
        index = 4
    else:
        index = 5
        
    if (maxFrac > args.frac):
        print "{}\t{}\t{:2.2f}\t{}\t{:2.2f}\t{}\t{}".format(vals[0], vals[1], maxFrac, index, readAlnLen, readLen, sinkLen)
    
    

    
