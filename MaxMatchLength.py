#!/usr/bin/env python

import Tools
import sys

samFile = open(sys.argv[1])
op = sys.argv[2]

totalLines = 0
for line in samFile:
    if (line[0] == '@'):
        continue
    else:
        sam = Tools.SAMEntry(line)
        maxMatch = 0
        if (sam.mapqv < 10):
            continue
        for i in range(len(sam.ops)):
            if (sam.ops[i] == op and sam.lengths[i] > maxMatch):
                maxMatch = sam.lengths[i]
        print maxMatch
    
        totalLines += 1

sys.stderr.write(str(totalLines) + "\n")
