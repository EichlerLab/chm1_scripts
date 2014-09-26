#!/usr/bin/env python

import AlignmentReaders
import sys

if (len(sys.argv) != 2):
    print "usage: IncompleteAlignments.py alignmetns.m1 . "
    sys.exit(1)
    
m4Name = sys.argv[1]

reader = AlignmentReaders.M4Reader(m4Name)

aln = reader.GetNext()
while (aln is not None):

    score = aln.qstart + (aln.qlen - aln.qend)
    if (aln.number == 0):
        print aln.qname + "\t" + str(score) + "\t" + str(aln.qstart) + "\t" + str(aln.qlen - aln.qend) +  "\t" + aln.ToString()
    aln = reader.GetNext()




