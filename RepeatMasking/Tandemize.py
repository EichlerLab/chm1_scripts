#!/usr/bin/env python

import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import argparse

ap = argparse.ArgumentParser(description="Make centromeric sinks from satellite monomers")
ap.add_argument("input", help="Input fasta")
ap.add_argument("output", help="Output fasta file.")
ap.add_argument("--length", help="Target length, repeat monomer until it is at least this long.", default=40000,type=int)


args = ap.parse_args()

inFasta = open(args.input)
outFasta = open(args.output, 'w')

for rec in SeqIO.parse(inFasta, "fasta"):
    recLen = len(rec.seq.tostring())
    nRec = (args.length / recLen ) 
    s=rec.seq.tostring() * (nRec + 1)
    expSeq = Seq.Seq(s)
    rec.seq = expSeq
    SeqIO.write(rec, outFasta, "fasta")

outFasta.close()

