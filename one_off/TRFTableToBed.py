#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])
for line in inFile.readlines():
    if (line[0] == '@'):
        titleVals = line[1:].split('/')
        if (len(titleVals) == 1):
            titleVals = [titleVals[0].strip(), 0, 0]
        else:
            titleVals[1] = int(titleVals[1])
            titleVals[2] = int(titleVals[2])
    else:
        vals = line.split()
        print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(titleVals[0], titleVals[1] + int(vals[0]), titleVals[1] + int(vals[1]), int(vals[1]) - int(vals[0]), vals[2], vals[3], vals[13])
        
    
