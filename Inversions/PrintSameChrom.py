#!/usr/bin/env python
import sys
if (len(sys.argv) < 2):
    print "usage: PrintSameChrom.py input.bed"
    sys.exit(1)
    
a = open(sys.argv[1])

prevVals = None
prevLine = ""
for line in a:
    vals = line.split()
    if (prevVals is not None and prevVals[0] == vals[0] and prevVals[3] == vals[3]):
        if (int(vals[1]) < int(prevVals[1])):
            sys.stdout.write(line)
            sys.stdout.write(prevLine)
        else:
            sys.stdout.write(prevLine)
            sys.stdout.write(line)
    prevLine = line
    prevVals = vals
