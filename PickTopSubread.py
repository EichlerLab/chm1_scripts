#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])

prn = ""
prl = 0
maxLen = 0
lineNumber = 0
for line in inFile:
    
    vals = line.split()
    i = vals[0].rfind("/")
    rn = vals[0][0:i]
    idx = vals[0][i+1:]
    (subStart, subEnd) = idx.split("_")
    subStart = int(subStart)
    subEnd   = int(subEnd)
    alnStart = int(vals[2])
    alnEnd   = int(vals[3])

    length = alnEnd - alnStart
    if (rn != prn):
        if (prn != ""):
            sys.stdout.write(mLine)
        mLine = line
        maxLen = length
    else:
        if (length > maxLen):
            mLine = line
            maxLength = length
    prn = rn        
    
