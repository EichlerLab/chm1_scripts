#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord


ap = argparse.ArgumentParser(description="Count masked characters.")
ap.add_argument("fastain", help="Input gap-bed file.")

args = ap.parse_args()


fastaFile = open(args.fastain)

for record in SeqIO.parse(fastaFile, "fasta") :


    seq = record.seq.tostring()
    nLower = seq.count("a") + seq.count("t") + seq.count("g") + seq.count("c")
    nN  = seq.count("N")
    if (nN == 0):
        print "{} {:2.2f}".format(record.id, float(nLower)/len(seq))
