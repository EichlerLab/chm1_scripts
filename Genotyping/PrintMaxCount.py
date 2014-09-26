#!/usr/bin/env python

import sys
import argparse
import Queue

ap = argparse.ArgumentParser()
ap.add_argument("infile", help="input")
ap.add_argument("--max", help="Keep this many.")

args= ap.parse_args()

inFile = open(args.infile)

prevRegion = ""
maxCount = 0
maxPat = ""
for line in inFile:
    vals = line.split()
    region = vals[0]
    if (region != prevRegion and prevRegion != ""):
        print prevRegion + "\t" + maxPat + "\t" + str(maxCount)
        maxCount = 0
    if (maxCount < int(vals[2])):
        maxCount = int(vals[2])
        maxPat = vals[1]
    prevRegion = vals[0]
print region + "\t" + maxPat + "\t" + str(maxCount)
