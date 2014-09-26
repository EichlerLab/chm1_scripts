#!/usr/bin/env python

import sys
if (len(sys.argv) < 3):
    print "usage: input.sam blacklist.txt [min]"
    sys.exit(0)
    
inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
minVal = 5
if (len(sys.argv) >= 4):
    minVal = int(sys.argv[3])
    
matches = {}
for line in inFile:
    if (line[0] != '@'):
        vals = line.split()
        if (vals[5] == "30M"):
            if (vals[9] not in matches):
                matches[vals[9]] = 0
            matches[vals[9]] += 1

for k in matches.keys():
    if (matches[k] >= minVal):
        outFile.write("{}\t{}\n".format(k, matches[k]))
        
    
