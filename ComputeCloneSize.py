#!/usr/bin/env python

import pysam
import numpy as np

import argparse
import sys

ap = argparse.ArgumentParser(description="Compute expected insert size")
ap.add_argument("bam", help="Input bam file.")
ap.add_argument("-n", help="Num samples", type=int, default=100000)
ap.add_argument("-M", help="Max isize", type=int, default=10000)
ap.add_argument("-m", help="Max isize", type=int, default=100)
ap.add_argument("--all", help="Print all to this file", default=None)
args = ap.parse_args()

bamFile = pysam.Samfile(args.bam, 'rb')

if (args.all is not None):
    allFile = open(args.all, 'w')

reads = {}
nUsed = 0
spans = []
nproc =0

for aln in bamFile.fetch():
    if (aln.qname not in reads):

        if (np.random.randint(0,50) < 1 and abs(aln.isize) < args.M and abs(aln.isize) > args.m ):
            spans.append(abs(aln.isize))
            reads[aln.qname]=True
    if (args.all is not None):
        allFile.write(str(aln.isize) + "\n")
        
    if (len(spans) >= args.n):
        break

    nproc += 1
    if (nproc % 10000 == 0):
        sys.stderr.write(str(nproc) + "\t" + str(len(spans)) + "\n")
        
npspans = np.asarray(spans)
npspans.sort()
print "Num: " + str(len(npspans))
print "Median: " + str(npspans[len(npspans)/2])
print "Mean: " + str(np.mean(npspans))
print "SD: " + str(np.std(npspans))
print "max: " + str(np.max(npspans))
print "max: " + str(np.max(spans))



