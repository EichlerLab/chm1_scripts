#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

for line in inFile:
    vals = line.split()
    if (vals[0] == "query"):
        continue
    gapRatio = float(vals[2]) / float(vals[1])
#    matchRatio = float(vals[3]) / float(vals[4])
    matchRatio = float(int(vals[3]) - 8000) / (float(vals[4])-8000)
    minRatio = float(int(vals[4]) - int(vals[2]))/int(vals[4])
    if (gapRatio < 0.25 and matchRatio > 0.9):
        line = line.strip()
        print line + "\t" + str(gapRatio) + "\t" +str(matchRatio)
