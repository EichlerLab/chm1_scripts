#!/usr/bin/env python

import sys
prevVals = [None, None, None]
totalCoverage = 0
nCoverage = 0
for line in sys.stdin.readlines():
	vals = line.split()
	if (vals[0:3] != prevVals[0:3] and prevVals[0] is not None):
		coverage = totalCoverage / nCoverage	
		prevVals[-1] = coverage 
		print "\t".join([str(i) for i in prevVals])
		totalCoverage = 0
		nCoverage = 0
	totalCoverage += float(vals[-1])
	nCoverage += 1.0
	prevVals = vals
