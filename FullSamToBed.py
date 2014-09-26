#!/usr/bin/env python

import Tools
import sys

samFile = open(sys.argv[1])

for line in samFile:
    if (line[0] == '@'):
        continue
    sam = Tools.SAMEntry(line)
    print sam.tName + "\t" + str(sam.tPos) + "\t" + str(sam.tEnd) + "\t" + line
