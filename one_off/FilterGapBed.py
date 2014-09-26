#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for line in inFile:
    vals = line.split()
    support = len(vals[7].split(';'))
    if (support > 2):
        outFile.write(line)
