#!/usr/bin/env python

import sys
if (sys.argv[1] == "-"):
    inFile = sys.stdin
else:
    inFile = open(sys.argv[1])

v=False
if (len(sys.argv) > 1 and sys.argv[2] == "-v"):
    v=True
    
for line in inFile:
    vals = line.split()
    reads = vals[3].split(';')
    readNames = {}
    for read in reads:
        base = "/".join(read.split("/")[0:2])
        if (base in readNames):
            readNames[base] = False
        else:
            readNames[base] = read
    sup = []
    for readName in readNames.values():
        if (readName != False):
            sup.append(readName)
    readNames = ";".join(sup)
    vals[3] = readNames
    vals[4] = str(len(sup))
    print "\t".join(vals)
    
