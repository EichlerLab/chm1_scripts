#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq

import sys

inFile = open(sys.argv[1])
for rec in SeqIO.parse(inFile, "fasta"):
    i = 0
    l = len(rec.seq)
    while (i < l):
        j=i
        while (j < l and rec.seq[j] == 'N'):
            j += 1
        if (j - i > 100000):
            print rec.id + "\t" + str(i) + "\t" + str(j)
        i = j+1
