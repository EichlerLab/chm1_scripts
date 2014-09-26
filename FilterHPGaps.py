#!/usr/bin/env python

import sys
bed = open(sys.argv[1])
nucs = ['A', 'C', 'G', 'T']
for line in bed:
    vals = line.split()
    seq = vals[5]
    counts = [float(seq.count(n))/len(seq) for n in nucs]
    printLine = True
    for c in counts:
        if (c >= 0.8):
            printLine = False
            break
    if (printLine):
        print line.strip()
            
