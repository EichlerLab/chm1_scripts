#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

for line in inFile:
    vals = line.split()
    if (vals[-1] == "left"):
        print vals[0] + "\t" + str(max(0,int(vals[1]))) + "\t" + str(int(vals[1]) + 300) + "\t" + vals[6]
    else:
        print vals[0] + "\t" + str(max(0,int(vals[2]) - 300)) + "\t" + vals[2] + "\t" + vals[6]
