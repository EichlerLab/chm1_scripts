#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("bedfile", help="File with allowable regions.")
ap.add_argument("countsFile", help="Original counts file.")
ap.add_argument("countsOut", help="Filtered counts file.")

args = ap.parse_args()

bedFile = open(args.bedfile)
queries = { "/".join(line.split()[0:3])  : True for line in bedFile }
countsFile = open(args.countsFile)
countsOut = open(args.countsOut, 'w')

for line in countsFile:
    vals = line.split()
    if (vals[0] in queries):
        countsOut.write(line)
