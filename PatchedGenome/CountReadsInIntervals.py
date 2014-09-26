#!/usr/bin/env python

import sys
import bisect

bedIn = open(sys.argv[1])
countFile = open(sys.argv[2])
bedOut = open(sys.argv[3], 'w')


#
# Build a lookup table for finding intervals.
#
start = {}
end   = {}
counts = {}
nIntervals = 0
for line in bedIn:
    vals = line.split()
    chrom = vals[0]
    if (chrom not in start):
        start[chrom] = []
        end[chrom] = []
        counts[chrom] = []        
    start[chrom].append(int(vals[1]))
    end[chrom].append(int(vals[2]))
    counts[chrom].append(0)
    nIntervals += 1
#sys.stderr.write("Done building start, {}\n".format(nIntervals))

nCounts = 0
nIncremented = 0

total = 0
for line in countFile:

    vals = line.split()
    if (len(vals) < 2):
        continue
    chrom = vals[0]
    pos = int(vals[1])
    count = int(vals[3])
    
    if (chrom in start):
        i = bisect.bisect_left(start[chrom], pos) - 1
#        if (i < len(start[chrom])-1):
#            pdb.set_trace()
#            print "before {}-{}, pos: {} after: {}-{}".format(start[chrom][i-1], end[chrom][i-1], pos, start[chrom][i+1], end[chrom][i+1])
        if (i >=0 and i < len(start[chrom]) and (start[chrom][i] <= pos and end[chrom][i] > pos) or (i < len(start[chrom])-1 and start[chrom][i+1] == pos)):
            counts[chrom][i] += count
            total += count
            nIncremented += 1
    nCounts += 1
#    if (nCounts % 100000 == 0):
#        sys.stderr.write("Incremented {} of {}\n".format(nIncremented, nCounts))
#sys.stderr.write("Searched {} counts\n".format(nCounts))
print str(total)
