#!/usr/bin/env python
import sys
if (len(sys.argv) != 2):
    print "usage: TopHitOnly.py in.sam out.sam"
    sys.exit(1)
inName = sys.argv[1]
outName = sys.argv[2]

inFile = open(inName, 'r')
outFile =open(outName, 'w')
prevAlign = ""
for line in inFile:
    if (len(line) > 0):
        if (line[0] == '@'):
            outFile.write(line)
        else:
            vals = line.split()
            if (vals[0] != prevAlign):
                outFile.write(line)
            prevAling = vals[0]
                
        
