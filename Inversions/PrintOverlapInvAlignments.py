#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="print overlapping alignments.")
ap.add_argument("input", help="input paired read bed file. adjacent lines are from the same read.")
ap.add_argument("--delta", help="min percent overlap, in range (0,1)", default=0.9, type=float)
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
    if (overlapPercent > args.delta):
        doPrint = True
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
            sys.stdout.write(l2)
