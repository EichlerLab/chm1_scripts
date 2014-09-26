#!/usr/bin/env python

import sys
if (len(sys.argv) < 2):
    print "Usage: print_bar_pos.py query_file "
    print " prints the indices of the boundaries between loci in the query file."
    sys.exit(1)
    
inFile = open(sys.argv[1])
pos = 0
vals = [ line.split() for line in inFile.readlines() ]
for i in range(0,len(vals)-1):
    if (vals[i][0] != vals[i+1][0] or i == len(vals)-2):
        print str(i+1)
