#!/usr/bin/env python

import sys
if (len(sys.argv) < 4):
    print "usage: RemoveBlacklist.py file.counts blacklist.txt file.new.counts"
    sys.exit(0)
    
countsFile = open(sys.argv[1])
blacklistFile = open(sys.argv[2])
countsOutFile = open(sys.argv[3], 'w')

blacklist = {}
for line in blacklistFile:
    vals = line.split()
    blacklist[vals[0]] = True

for line in countsFile:
    vals = line.split()
    if (vals[1] not in blacklist):
        countsOutFile.write(line)
        
