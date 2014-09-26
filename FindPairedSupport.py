#!/usr/bin/env python


import sys
if (len(sys.argv) < 3):
    print "usage: FindPairedSupport.py table left-isect.bed right-isect.bed\n"
    sys.exit(0)
table = open(sys.argv[1])
left  = open(sys.argv[2])
right = open(sys.argv[3])

leftDict = {}
for line in left:
    vals = line.split()
    title = vals[0] + ":" + vals[1] + "-" + vals[2]
    leftDict[title] = int(vals[-1])
    
rightDict = {}
for line in right:
    vals = line.split()
    title = vals[0] + ":" + vals[1] + "-" + vals[2]
    rightDict[title] = int(vals[-1])

for line in table:
    vals = line.split()
    leftTitle = vals[0] + ":" + vals[1] + "-" + vals[2]
    rightTitle = vals[3] + ":" + vals[4] + "-" + vals[5]
    if (leftTitle in leftDict and rightTitle in rightDict):
        line = line.strip()
        print line + "\t" + leftTitle + "\t" + rightTitle
