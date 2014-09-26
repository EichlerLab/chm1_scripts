#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio import Seq

inFile = open(sys.argv[1])

for rec in SeqIO.parse(inFile, "fasta"):
    i = 0
    l = len(rec.seq)
    while (i < l):
        j = i;
        while (rec.seq[j] == 'N'):
            j+=1

        i = j+1
        while (i+1 < l and rec.seq[i+1] != 'N'):
            i+=1
