#!/usr/bin/env python

from Bio import SeqIO
import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for seq in SeqIO.parse(inFile, "fasta"):
    SeqIO.write(seq, outFile, "fasta")

outFile.close()
