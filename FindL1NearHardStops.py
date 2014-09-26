#!/usr/bin/env python

import sys
import argparse
import Tools
import os

ap = argparse.ArgumentParser(description="Find L1s near alignments.")
ap.add_argument("l1", help="L1 alignments in m1 format.")
ap.add_argument("hs", help="Input hard-stop cluster bed file.  The 4th field is the reads that align")

args = ap.parse_args()


l1File = open(args.l1)
l1Alignments = {}
    
for line in l1File:
    vals = line.split()
    if (int(vals[3]) == 1):
        qlen = int(vals[11])
        qs   = qlen - int(vals[10])
        qe   = qlen - int(vals[9])
        vals[9] = qs
        vals[10] = qe
    l1Alignments[vals[0]] = vals[1:]

hsFile = open(args.hs)
for line in hsFile:
    vals = line.split()
    reads = vals[3].split(";")
    for read in reads:
        if (read in l1Alignments):
            print "\t".join(vals[0:3]) + "\t" + read
