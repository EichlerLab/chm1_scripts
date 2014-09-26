#!/usr/bin/env python

import sys

inFile = open(sys.argv[1])

while True:
    l1 = inFile.readline()
    l2 = inFile.readline()
    if (l1 == "" or l2 == ""):
        break
    v1 = l1.split()
    v2 = l2.split()
    if (v1[8] != v2[8]):
        sys.stdout.write(l1)
        sys.stdout.write(l2)
        
    
