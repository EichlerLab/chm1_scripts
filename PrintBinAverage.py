#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="print average coverage in a region.")
ap.add_argument("--inFile", help="Read input from.", default="stdin")
ap.add_argument("--outFile", help="Write here.", default="stdout")
ap.add_argument("--window", help="Window size", type=int,default=50)
args = ap.parse_args()


if (args.inFile == "stdin"):
    inFile = sys.stdin
else:
    inFile = open(args.inFile)

prevLine = ""
prevChrom = ""
span = 0
binCoverage = 0.0
prevStart = None
#while (line = inFile.readline()):
start = 0
while (inFile):
    line = inFile.readline()
    
    if (line != ""):
        vals = line.split()
        if (vals[0] == prevChrom):
            binCoverage += int(vals[2])
            span += 1
    else:
        break
    # print on chromosome boundaries or window boundaries or eof
    if ((prevChrom != vals[0] and prevChrom != "") or span == args.window or line == ""):
        if (span != 0):
            averageCoverage = binCoverage / span
            print "{}\t{}\t{}\t{:.2f}".format(prevChrom, start, start + span, averageCoverage)

        span = 1
        binCoverage = float(int(vals[2]))
        start = int(vals[1])
  
    prevChrom = vals[0]
    
    
    
        
