#!/usr/bin/env python

import sys

fh = open(sys.argv[1])
for line in fh:
    vals = line.split()
    leftReads = vals[3].split(';')
    leftSupport = {}
    rightReads = vals[8].split(';')
    
    for leftRead in leftReads:
        base = "/".join(leftRead.split("/")[0:2])
        leftSupport[base] = True

    rightSupport = {}
    for rightRead in rightReads:
        base = "/".join(rightRead.split("/")[0:2])
        rightSupport[base] = True
        
    if (len(leftSupport) > 2 and len(leftSupport) < 20 and len(rightSupport) > 2 and len(rightSupport) < 20):
        sys.stdout.write(str(len(leftSupport)) + "\t" + line)
        
