#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord


ap = argparse.ArgumentParser(description="Print sequences that are not well masked to an output file.")
ap.add_argument("bedin", help="Input gap-bed file.")
ap.add_argument("bedpass", help="Output gap-bed file, passing filter.")
ap.add_argument("bedfail", help="Output sequences with too low masked content to this file.")
ap.add_argument("--minMasked", help="Minimum masking content", type=float, default=0.80)

args = ap.parse_args()


bedIn = open(args.bedin)
bedPass = open(args.bedpass, 'w')
bedFail = open(args.bedfail, 'w')

for line in bedIn:
    vals = line.split()
    seq = vals[5]
#    nUpper = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
    nLower = seq.count("a") + seq.count("t") + seq.count("g") + seq.count("c")

    ratio = (float(nLower) / len(seq))
    if (ratio >= args.minMasked):
        bedPass.write(line)
    else:
        bedFail.write(line)
    
