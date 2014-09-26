#!/usr/bin/env python

import sys

def Close(a,b,l):
    return (abs(a-b) < l)

inFile = open(sys.argv[1])
leftFile = open(sys.argv[2], 'w')
rightFile = open(sys.argv[3], 'w')

while True:
    l1 = inFile.readline()
    l2 = inFile.readline()

    if (l1 == "" or l2 == ""):
        break
    v1 = l1.split()
    v2 = l2.split()
    i1 = [v1[0], int(v1[1]), int(v1[2])]
    i2 = [v2[0], int(v2[1]), int(v2[2])]

    # check for close breakpoints.
    if (Close(i1[1], i2[1], 100) and not Close(i1[2], i2[2], 100) and v1[8] != v2[8]):
        leftFile.write(l1)
        leftFile.write(l2)        
        
    if (Close(i1[2], i2[2], 100) and not Close(i1[1], i2[1], 100) and v1[8] != v2[8]):
        rightFile.write(l1)
        rightFile.write(l2)        
        
        
    
