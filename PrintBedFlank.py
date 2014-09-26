#!/usr/bin/env python

import sys
import argparse
import Tools
ap = argparse.ArgumentParser(description="Print a window of either the left side or right side of a bed interval.")
ap.add_argument("--input", help="input file, or stdin", default="stdin")
ap.add_argument("--output", help="output file, or stdout", default="stdout")
#ap.add_argument("--left", help="Print window on left side.", action='store_true', default=False)
#ap.add_argument("--right", help="Print window on right side.", action='store_true', default=False)
ap.add_argument("--fai", help="Read this fasta index to double check coordinates with chromosome size.", default=None)
ap.add_argument("--internal", help="Use this much on side of the interval.", default=100, type=int)
ap.add_argument("--external", help="Add this much on outside of interval, for example, subtract from left interval.", default=0, type=int)

args = ap.parse_args()
if (args.input == "stdin"):
    inFile = sys.stdin
else:
    inFile = open(args.input)

if (args.output == "stdout"):
    outFile = sys.stdout
else:
    outFile = open(args.output, "w")

if (args.fai is not None):
    fai = Tools.ReadFAIFile(args.fai)


def BoundedIncrease(val, i, chrom, fai):
    if (fai is None or chrom not in fai):
        return val + i
    else:
        return min(val+i,fai[chrom][0])

numReordered = 0    
for line in inFile.readlines():
    vals = line.split()
    for i in (0,8):
        start = int(vals[i+1])
        end   = int(vals[i+2])

        if (vals[i+7] == "left"):
            vals[i+1] = str(max(0, start - args.external))
            vals[i+2] = str(BoundedIncrease(start, args.internal, vals[0], args.fai))
        if (vals[i+7] == "right"):
            vals[i+1] = str(max(0, start - args.internal))
            vals[i+2] = str(BoundedIncrease(end, args.external, vals[0], args.fai))
    if (vals[0] == vals[8] and int(vals[1]) >= int(vals[9])):
        vals2 = vals[8:16] + vals[0:8]
        if (len(vals2)  != len(vals)):
            print str(vals)
            print str(vals2)
            sys.exit(1)
        vals = vals2
        numReordered += 1
    outFile.write("\t".join(vals) + "\n")
sys.stderr.write("Reordered " + str(numReordered) + "\n")
if (outFile is not sys.stdout):
    outFile.close()
        
        

        
    
    
