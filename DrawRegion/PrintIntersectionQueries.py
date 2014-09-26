#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])

for line in inFile:
    vals = line.split()
    intv = vals[0].split("/")[3].split("_")
    start = int(intv[0])
    end   = int(intv[1])
    print vals[0] + "\t4000\t"+ str(4000+(end-8000))
