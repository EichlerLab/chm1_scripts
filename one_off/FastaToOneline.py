#!/usr/bin/env python

from Bio import SeqIO

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for seq in SeqIO.parse(inFile, "fasta"):
    outFile.write(seq.id+"\n")
    outFile.write(str(seq.seq) +"\n")
    
