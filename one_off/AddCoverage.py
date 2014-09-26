#!/usr/bin/env python

import pysam
import sys


inFileName = sys.argv[1]
bamFile = pysam.Samfile(sys.argv[2])
outFile = open(sys.argv[3], 'w')

if (inFileName == "stdin"):
    inFile = sys.stdin
else:
    inFile = open(inFileName)

for line in inFile.readlines():
    region = line.split()
    sup = bamFile.count(region[0], (int(region[2]) + int(region[1])) / 2 - 1, (int(region[2]) + int(region[1])) / 2)
    outFile.write(line.strip() + "\t" + str(sup) + "\n")
outFile.close()
