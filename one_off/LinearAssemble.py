#!/usr/bin/env python

import sys
inFileName = sys.argv[1]


kmerTableFile = open(inFileName, 'r')

curSeq = ""
seqIndex = 1
prevTitle = ""
for line in kmerTableFile:
    vals = line.split()
    vals[1] = vals[1].upper()
    if (curSeq == ""):
        curSeq = vals[1]
    else:
        l = len(vals[1])
        suf = curSeq[-l+1:]
        pre = vals[1][:-1]

        if (suf == pre):
            prevTitle = vals[0]            
            curSeq += vals[1][-1]
        else:
            print ">" + prevTitle
            print curSeq
            curSeq = vals[1]
            seqIndex += 1
