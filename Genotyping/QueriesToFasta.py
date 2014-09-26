#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for line in inFile:
    vals = line.split()
    outFile.write(">"+vals[0]+"\n")
    outFile.write(vals[1]+"\n")
outFile.close()
