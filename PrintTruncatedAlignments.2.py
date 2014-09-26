#!/usr/bin/env python

import Tools
import sys
import argparse
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq

ap = argparse.ArgumentParser(description="Print alignments (much) shorter than the expected read length.")
ap.add_argument("sam", help="Input SAM file (s)", nargs="+")
ap.add_argument("--minTrunc", help="Minimal truncated length.", type=int, default=500)
ap.add_argument("--minAlign", help="Minimal aligned length.", type=int, default=1000)
ap.add_argument("--out", help="Write output here.", default=None)

args = ap.parse_args()


starti = 5
endi   = 6
flagi  = 1
tposi  = 3
qnamei = 0
tnamei = 2
seqi   = 8
prevSamEntry = None
prevWasShort = False

if (args.out is not None):
    out = open(args.out, 'w')
else:
    out = sys.stdout


for samFileName in args.sam:
    samFile = open(samFileName, 'r')
    sys.stderr.write(samFileName + "\n")
    for line in samFile:
        if (len(line) > 0 and line[0] == '@'):
            continue
        aln = Tools.SAMEntry(line)
        # only process primary alignments
        if (aln.flag & 256 != 0):
            continue
        if (aln.qEnd - aln.qStart < args.minAlign):
            continue
        
        (barcode, zmw, coords) = Tools.ParseReadTitle(aln.title)
        delta = (coords[1] - coords[0]) - (aln.qEnd - aln.qStart )

        if (coords[1] - aln.qEnd > aln.qStart - coords[0]):
            side = "right"
        else:
            side = "left"
        if ( abs(delta) > args.minTrunc):
            out.write( aln.tName + "\t" + str(aln.tStart) + "\t" + str(aln.tEnd) + "\t" + str(delta) + "\t" + aln.title + "\t" + str(delta) + "\t" + side + "\n")
if (out != sys.stdout):
    out.close()
