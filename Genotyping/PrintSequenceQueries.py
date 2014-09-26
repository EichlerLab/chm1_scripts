#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

for line in inFile:
    vals = line.split()
    name = "/".join(vals[0:3])
    start = 4000
    end = start + int(vals[4])
    print name + "\t" + str(start) + "\t" + str(end)
