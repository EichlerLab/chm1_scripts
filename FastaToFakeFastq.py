#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

import sys
inFile = open(sys.argv[1],'r')
outFile = open(sys.argv[2], 'w')

for rec in SeqIO.parse(inFile, "fasta"):
    rec.letter_annotations['phred_quality'] =  [ 20 ]*len(rec.seq)
    SeqIO.write(rec, outFile, "fastq-sanger")

outFile.close()
