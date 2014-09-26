#!/usr/bin/env python

import sys
import Tools
import argparse

ap = argparse.ArgumentParser(description="print either stretched or not stretched alignments.")
ap.add_argument("sam", help="Input file.")
ap.add_argument("-v", help="Not stretched.", default=False, action='store_true')
ap.add_argument("--rate", help="Min stretched (0.05)", type=float, default=0.05)
ap.add_argument("--readlen", help="Use read length, not aligned length.", action='store_true', default=False)
args = ap.parse_args()
inFile = open(args.sam)
for line in inFile:
    if (line[0] == '@'):
        continue
    vals = line.split()
    if (len(vals) == 0):
        continue
    xs = int(Tools.GetKV("XS:i:", vals))
    xe = int(Tools.GetKV("XE:i:", vals))
    
    (readStart, readEnd) = vals[0].split("/")[1].split("_")
    readStart = int(readStart)
    readEnd   = int(readEnd)

    readLen = readEnd - readStart

    alnLen = xe - xs
    targetLen = int(vals[8])
    l = alnLen
    if (args.readlen):
        l = readLen

    targetRatio = float(targetLen) / float(l)
    op = None
    if (targetRatio > (1 + args.rate)):
        op = "del"
    if (targetRatio < (1/(1+args.rate))):
        op = "ins"
    if (op is not None and args.v == False):
        print vals[0] + "\t" + op + "\t{:2.2f}".format(targetRatio) 

    if (op is None and args.v == True):
        print vals[0]
        
    
    
