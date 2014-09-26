#!/usr/bin/env python

import argparse
import sys

actionSide = 1

ap = argparse.ArgumentParser(description="Modify a bed file.")
subparsers = ap.add_subparsers(help="commands")
parserSide = subparsers.add_parser("side", help="Print left or right side of bed.")
parserSide.add_argument("bed", help="Input bed file.")
parserSide.add_argument("--out", help="Output file (default is stdout)", default="/dev/stdout")
parserSide.add_argument("--delta", help="Print left + delta", default=10,type=int)
parserSide.add_argument("--right", help="Print right side, default is left", action='store_true', default=False)
parserSide.set_defaults(action=actionSide)

args = ap.parse_args()


bedFile = open(args.bed)
outFile = open(args.out, 'w')

if (args.action == actionSide):
    for line in bedFile:
        vals = line.split()
        if (args.right == True):
            vals[1] = str(max(int(vals[1]), int(vals[2]) - args.delta))
        else:
            vals[2] = str(min(int(vals[1]) + args.delta, int(vals[2])))
        outFile.write("\t".join(vals) + "\n")
outFile.close()




