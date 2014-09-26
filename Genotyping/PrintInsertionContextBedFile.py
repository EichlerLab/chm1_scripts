#!/usr/bin/env python

import sys
from Bio import SeqIO

if (len(sys.argv) <= 1):
    print "Usage: PrintInsertionContextBedFile.py file.fasta window. "
    print "This will output a bed file for the coordinates of the inserted sequence (not flanking)"
    sys.exit(1)

fastaFile = open(sys.argv[1])
window = int(sys.argv[2])



for seq in SeqIO.parse(fastaFile, "fasta"):
    print seq.id + "\t" + str(window) + "\t" + str(len(seq.seq)-window)
fastaFile.close()

