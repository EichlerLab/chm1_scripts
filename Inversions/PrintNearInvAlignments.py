#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="print overlapping alignments.")
ap.add_argument("input", help="input paired read bed file. adjacent lines are from the same read.")
ap.add_argument("--window", help="Maximum allowed window", default=1000, type=int)
ap.add_argument("--unique_subreads", help="Only print support from one subread.", action='store_true', default=False)
ap.add_argument("--opposite", help="Print if on opposite strand.", action='store_true', default=False)

args = ap.parse_args()

inFile = open(args.input)

reads={}

while True:
    l1 = inFile.readline()
    l2 = inFile.readline()
    if (l1 == "" or l2 == ""):
        break
    v1 = l1.split()
    v2 = l2.split()
    intv1 = [v1[0], int(v1[1]), int(v1[2])]
    intv2 = [v2[0], int(v2[1]), int(v2[2])]

    
    overlapStart = max(intv1[1], intv2[1])
    overlapEnd   = min(intv1[2], intv2[2])
    overlap=max(0, overlapEnd - overlapStart)
    overlapPercent = float(overlap) / (intv2[2] - intv2[1])

    subreadCoordinates = v1[3].split("/")[-1].split("_")
    subreadLength = int(subreadCoordinates[1]) - int(subreadCoordinates[0])

    v1ts = int(v1[1])
    v1te = int(v1[2])
    v2ts = int(v2[1])
    v2te = int(v2[2])
    
    v1qs = int(v1[6])
    v1qe = subreadLength - int(v1[7])
    v2qs = int(v2[6])
    v2qe = subreadLength - int(v2[7])


    if (v2te >= v1te):
        gap = min(v2ts - v1te, v2qs - v1qe)
    else:
        gap = min(v1qs - v2qe, v1ts - v2te)

    doPrint = True
    if (gap > args.window):
        doPrint = False
        continue

    if (args.unique_subreads):
        read = "/".join(v1[3].split('/')[0:2])
        if read in reads:
            doPrint = False
        else:
            reads[read] = True

    if (args.opposite):
        if (v1[8] == v2[8]):
            doPrint = False
                
    if (doPrint):
#        sys.stdout.write(l1)
        sys.stdout.write(l2)
