#!/usr/bin/env python

from Bio import SeqIO
import sys

inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')


for record in SeqIO.parse(inFile, "fasta"):
    nN = record.seq.count("N")
    if (nN > 0):
        continue
    else:
        outFile.write(record.id + "\n")
        
