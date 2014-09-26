#!/usr/bin/env python

import pysam
import argparse
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser(description="Print the nucleotide frequency of a pileup.")
ap.add_argument("--bam", help="Input bam file.", required=True)
ap.add_argument("--freq", help="Print the frequency to this file", required = True)
ap.add_argument("--plot", help="Plot a SVG of the frequency")
ap.add_argument("--max", help="Limit alignments processed", default=None, type=int)


args = ap.parse_args()


bamFile = pysam.Samfile(args.bam, 'rb')

m = [4]*256
m[ord('A')] = 0
m[ord('C')] = 1
m[ord('G')] = 2
m[ord('T')] = 3


idx = 0
sequences = {}
for seq in bamFile.header["SQ"]:
    print "adding ref of len " + str(seq["LN"])
    sequences[idx] = np.zeros((4,seq["LN"]))
    idx+=1

index = 0
for aln in bamFile.fetch():
    mat = sequences[aln.tid]
    tPos = aln.pos
    qPos = aln.qstart
    qseq = aln.query
    for op in aln.cigar:
        if (op[0] != 0 and op[0] != 1 and op[0] != 2):
            continue
        if (op[0] == 0):
            for i in range(op[1]):
                mat[m[ord(qseq[qPos])]][tPos] += 1
                qPos +=1
                tPos +=1
        elif (op[0] == 1):
            qPos += op[1]
        elif (op[0] == 2):
            tPos += op[1]
    index += 1
    if (index % 1000 == 0):
        print "processed: " + str(index)
    if (args.max is not None and args.max == index):
        break

outFile = open(args.freq, 'wb')
for seq in sequences:
    np.save(outFile, sequences[seq])

outFile.close()

