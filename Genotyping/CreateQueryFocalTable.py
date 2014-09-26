#!/usr/bin/env python

import sys
from Bio import SeqIO

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')



for seq in SeqIO.parse(inFile, "fasta") :
    name = seq.id
    vals = seq.id.split("/")
    start = int(vals[1])
    end   = int(vals[2])
    totalLength = len(seq.seq)
    prefix = totalLength - (end - start)
    focusStart = prefix
    focusEnd  = prefix + (end - start)
    outFile.write( name + "\t" + str(focusStart) + "\t" + str(focusEnd) + "\n")
