#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
names = {}
for line in inFile.readlines():
    vals = line.split()
    name = "|".join(vals[0:3])
    if (name in names):
        print "removing " + str(int(vals[2]) - int(vals[1])) + "\t" + "\t".join(vals[0:3])
        continue
    else:
        names[name] = True
        outFile.write(line)
