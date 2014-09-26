#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq


import sys

inFile = open(sys.argv[1])
seqs = {}

for rec in SeqIO.parse(inFile, "fasta"):
    fullName = rec.id
    baseName = "/".join(fullName.split("/")[0:2])
    l = len(rec.seq)
    if (baseName not in seqs):
        seqs[baseName] = l
    else:
        if (l > seqs[baseName]):
            seqs[baseName] = l

for k,v in seqs.iteritems():
    print str(k) + "\t" + str(v)
