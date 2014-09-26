#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

prevChrom = None
prevPos   = None
total = 0
nPos = 0
last = False
import pdb


for line in inFile:

    vals = line.split()
    if (len(vals) != 3):
        last=True
        endPos = prevPos
    else:
        chrom = vals[0]
        pos   = int(vals[1])
        cov   = int(vals[2])
    # detect when to print out an interval
    atBoundary = False
    if (vals[0] != prevChrom):
  	atBoundary = True
    if (pos % 500 == 0):
	atBoundary = True
    if (prevPos is not None and pos != prevPos + 1):
	atBoundary = True
    if (last == True):
	atBoundary = True 
    if (atBoundary and prevChrom is not None):
	endPos = prevPos
	if (prevPos != boundaryStart):
	    printStart = (boundaryStart/ 500)*500
	    printEnd  = printStart + 500
            print vals[0] + ":" + str(printStart) + "-" + str(printEnd) + "\t{:2.2f}".format(float(total)/(endPos - boundaryStart + 1))
        total = 0
    if (atBoundary or prevChrom is None):
	boundaryStart = pos

    prevPos = pos
    prevChrom = chrom
        
    total += cov    

    if (last):
        break
