#!/usr/bin/env python

import sys
import argparse
ap= argparse.ArgumentParser(description="Print various types of relations between pairs of alignments.")
ap.add_argument("input", help="hardtop output file, sorted by read name (field 4)")
ap.add_argument("--relation", help="Type of overlap allowed: cover|exclude|close ", default="exclude", choices=["cover","exclude", "close"])
ap.add_argument("--delta", help="Max difference between breakpoints.", default=50000)
ap.add_argument("--parity", help="same, or opposite orientation", choices=["same", "opposite", "any"], default="any")

args = ap.parse_args()

inFile = open(args.input)

prevVals = None
prevLine = ""
for line in inFile:
    vals = line.split()
    if (prevVals is not None and prevVals[3] == vals[3]):
        readIntv = [int(i) for i in vals[3].split("/")[2].split("_")]
        readLen  = readIntv[1] - readIntv[0]

        (aStart, aEnd) = (int(vals[6]), readLen - int(vals[7]))
        (bStart, bEnd) = (int(prevVals[6]), readLen - int(prevVals[7]))
        aInsideB = (aStart >= bStart and aEnd <= bEnd)
        bInsideA = (bStart >= aStart and bEnd <= aEnd)
        if (len(vals) == 8):
            parity = "any"
        elif (vals[8] == prevVals[8]):
            parity = "same"
        else:
            parity = "opposite"
            
        if (args.relation == "exclude"):
            if (aInsideB == False and bInsideA == False):
                if (args.parity == "any" or parity == args.parity):
                    sys.stdout.write(prevLine)
                    sys.stdout.write(line)
        elif (args.relation == "cover"):
            if (aInsideB == True or bInsideA == True):
                if (args.parity == "any" or parity == args.parity):
                    sys.stdout.write(prevLine)
                    sys.stdout.write(line)
        elif (args.relation == "close"):
            if (min(abs(aStart - bStart), abs(aStart - bEnd), abs(aEnd - bStart), abs(aEnd - bEnd)) < args.delta):
                if (args.parity == "any" or parity == args.parity):
                    sys.stdout.write(prevLine)
                    sys.stdout.write(line)
                
                
    prevVals = vals
    prevLine = line
