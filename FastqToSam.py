#!/usr/bin/env python

import sys

index = 0
while (sys.stdin):
    name = sys.stdin.readline()
    if (name == ""):
        break
    seq  = sys.stdin.readline()
    sep  = sys.stdin.readline()
    qual = sys.stdin.readline()
    tlen = len(seq)
    sys.stdout.write(name[1:].strip() + "\t4\t*\t0\t0\t*\t*\t0\t" + str(tlen) + "\t"+seq.strip() + "\t" + qual)
    if (index % 100000 == 0):
        sys.stderr.write(str(index) + "\n")
    index +=1
    
    
