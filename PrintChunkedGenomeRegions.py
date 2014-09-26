#!/usr/bin/env python

import argparse
import sys
import Tools

ap = argparse.ArgumentParser(description="print reigons using a .fai file")
ap.add_argument("fai", help="fai file")
ap.add_argument("--chrs", help="Print for these chromosomes,default is all.", nargs="+", default=None)
ap.add_argument("--start", help="Start position", default=0, type=int)
ap.add_argument("--end", help="End pos",default=0, type=int)
ap.add_argument("--size", help="Print this size of chunk.", default=10000, type=int)
ap.add_argument("--stride", help="Skip this many bases between chunks.", default=5000,type=int)

args = ap.parse_args()

fai = Tools.ReadFAIFile(args.fai)
if (args.chrs is None):
    chrs = fai.keys()

else:
    chrs = args.chrs
rangeStrs = []
for chrom in chrs:
    chrLen = fai[chrom][0]
    start = args.start
    if args.end == 0:
        end = chrLen
    else:
        end = args.end
        
    for i in range(0,chrLen, args.stride):
        print chrom + ":" + str(i) + "-" + str(min(chrLen,i+args.size))




