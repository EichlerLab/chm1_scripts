#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for line in inFile:
    vals = line.split()
    for i in range(3,len(vals)):
        if (vals[i][0:3] == "chr"):
            break
    
    if (i < len(vals)):    
        titleVal = vals[i]
        if (titleVal.find(vals[0]) == -1):
            print "MISSING " + line
        else:
            outFile.write(line)
outFile.close()

