#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

for line in inFile:
    vals = line.split()
    seq = vals[5].upper()
    l = len(seq)
    nucs = ['A', 'T', 'C', 'G']
    fracs = [float(seq.count(i))/l for i in nucs]
    seqBias = ""
    for i in range(0,len(nucs)):
        if (fracs[i] > 0.35):
            seqBias += nucs[i]
    if (seqBias == ""):
        seqBias = "N"
    sys.stdout.write("\t".join(vals[0:3]) + "\t" + seqBias + "\t" + vals[4]   + "\n")
    
        
        
    
