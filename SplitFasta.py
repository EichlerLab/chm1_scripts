#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Split fasta files in different ways.")

ap.add_argument("--input", help="Input file (stdin)", default="/dev/stdin")
ap.add_argument("--output", help="Output file base", default="split")
ap.add_argument("--size", help="Maximum size of each smaller file.", default=None, type=int)
ap.add_argument("--nReads", help="Maximum number of sequences per file.", default=None, type=int)

args = ap.parse_args()

inFile = open(args.input)

nReads = 0
nBases = 0


curRead = 0
curSize = 0


outFile = None
fileIndex = 0

for seq in SeqIO.parse(inFile, "fasta"):
    createNewFile = False
    if (outFile is None):
        createNewFile = True
    if (args.size is not None and curSize + len(seq) > args.size):
        createNewFile = True
    if (args.nReads is not None and curRead == args.maxN):
        createNewFile = True

    if (createNewFile):
        if (outFile is not None):
            outFile.close()
        fileIndex += 1
        sys.stderr.write("{} {} {}\n".format(fileIndex, curRead, curSize))
        outFile = open(args.output + "." + str(fileIndex) + ".fasta", 'w')
        curRead = 0
        curSize = 0
    SeqIO.write(seq, outFile, "fasta")
    curRead += 1
    curSize += len(seq)
    
inFile.close()
