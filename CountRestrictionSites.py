#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse

ap = argparse.ArgumentParser(description="Count some statistics about restriction cut sites from either reads or contigs.")
ap.add_argument("sequence", help="Input sequence")
ap.add_argument("--motifs", help="Restriction recognitition sequences.", nargs='+')
ap.add_argument("--printAll", help="Print all counts", action='store_true', default=False)


opts = ap.parse_args()

h = open(opts.sequence, 'r')
counti = { m: 0 for m in opts.motifs }
existi = { m: 0 for m in opts.motifs }
nseq = 0
for seq in SeqIO.parse(h, "fasta"):
    mCount = {}
    for motif in opts.motifs:
        cm = seq.seq.upper().count(motif)
        mCount[motif] = cm
        counti[motif] += cm
        if (cm != 0):
            existi[motif] += 1
        
    if (opts.printAll):
        print mCount
    nseq += 1

print "nseq " + str(nseq)
sys.stdout.write("total motifs: ")
for c in counti:
    sys.stdout.write(" " + str(counti[c]))
sys.stdout.write("\n")

sys.stdout.write("with motif: ")
for c in existi:
    sys.stdout.write(" " + str(existi[c]))
sys.stdout.write("\n")
