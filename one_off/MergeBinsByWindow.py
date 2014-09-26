#!/usr/bin/env python
import argparse
ap = argparse.ArgumentParser(description="Merge bins in the same window.")
ap.add_argument("bed", help="Input bed")
ap.add_argument("out", help="Output bed")
ap.add_argument("--window", help="Window size", type=int, default=5000)

args = ap.parse_args()

bedIn = open(args.bed)
bedOut = open(args.out, 'w')

curChr = ""

for line in bedIn:
    vals = line.split()
    if (vals[0] != curChr):
        binStartChr = vals[0]
        binStartPos = int(vals[1])
        binEndPos   = int(vals[2])
    else:
        s = int(vals[1])
        e = int(vals[2])
        if (e - binStartPos < args.window):
            binEndPos = e
        else:
            bedOut.write("{}\t{}\t{}\n".format(binStartChr, binStartPos, binEndPos))
            binStartPos = s
            binEndPos   = e
            binStartChr = vals[0]
    curChr = vals[0]
                         
