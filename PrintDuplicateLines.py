#!/usr/bin/env python
import sys
import argparse
ap = argparse.ArgumentParser(description="Print duplicate hard-stop lines after they have been sorted by read name.")
ap.add_argument("input", help="Input lines")
ap.add_argument("out", help="Output lines")
ap.add_argument("--orientation", help="Orientation to print, ll, lr, rr, any, same", default="any")
ap.add_argument("-v", help="Opposite, print unique lines.", action='store_true', default=False)

args = ap.parse_args()
inFile = open(args.input, 'r')
outFile = open(args.out, 'w')

prevRead = ""
prevSide = ""
prevLine = ""
prevChrom = ""
for line in inFile:
    vals = line.split()
    curSide = vals[-1]
    
    if (vals[3] == prevRead):
        if (args.orientation == "any" or
            (args.orientation == "ll" and curSide == "left" and prevSide == "left") or
            (args.orientation == "rr" and curSide == "right" and prevSide == "right") or
            (args.orientation == "lr" and curSide != prevSide ) or
            (args.orientation == "same" and curSide == prevSide )):
            if (args.v == False and vals[0] == prevChrom):
                outFile.write(prevLine)
                outFile.write( line)
                
        prevRead = ""
        prevChrom = ""
    else:
        if (args.v == True and prevLine != ""):
            outFile.write(prevLine)
        prevRead = vals[3]
        prevLine = line
        prevSide = curSide
        prevChrom = vals[0]
        
outFile.close()
