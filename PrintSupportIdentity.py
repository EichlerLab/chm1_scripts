#!/usr/bin/env python

import Align
import sys
import pdb

inFile = open(sys.argv[1])
nMatch = 0
nTest = 0
for line in inFile.readlines():
    vals = line.split()
    seqs = vals[6].split(';')
    if (len(seqs) == 2):
        (qopt, topt, score) = Align.SWAlign(seqs[0], seqs[1], indel=0,mismatch=0)
        ident = float(score) / max(len(seqs[0]), len(seqs[1]))
        if (ident > 0.70):
            print line
            nMatch += 1
        nTest +=1
        if (nTest % 1000 == 0):
            sys.stderr.write(str((nTest, nMatch)) + "\n")
    else:
        print line
#        print str((seqs[0], seqs[1], score))
#        print "{} {} {}".format(
