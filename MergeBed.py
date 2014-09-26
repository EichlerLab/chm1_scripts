#!/usr/bin/env python

import sys
bedIn = open(sys.argv[1], 'r')
bedOut = open(sys.argv[2], 'w')

pc = None
for line in bedIn:
    vals = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end   = int(vals[2])
    if (pc != None and
        chrom == pc[0]
        and ((pc[0] <= start and pc[1] > start) or
             (pc[1] >= start and pc[1] > end) or
             (pc[0] <= start and pc[1] > end))):
        pass
    else:
        bedOut.write(line)
    pc = (chrom, start, end)
