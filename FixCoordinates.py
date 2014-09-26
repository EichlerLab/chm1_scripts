#!/usr/bin/env python

import Tools
import argparse
import sys
ap = argparse.ArgumentParser(description="Fix bed coordinates that may have gone past the end of a reference contig to be flush with the end.")
ap.add_argument("input", help="Input BED file.")
ap.add_argument("output", help="Output BED file.")
ap.add_argument("fai", help="Reference fai file.")
ap.add_argument("--col", help="Column to repair (2)", type=int, default=2)
args = ap.parse_args()

bedIn = open(args.input)
bedOut = open(args.output, 'w')

fai = Tools.ReadFAIFile(args.fai)

for line in bedIn:
    vals = line.split()
    if (vals[0] in fai):
	vals[args.col] = str(min(int(vals[args.col]), fai[vals[0]][0]))
    bedOut.write('\t'.join(vals)+'\n')

    
        
        
