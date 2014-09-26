#!/usr/bin/env python

import Tools
import argparse
import sys

ap = argparse.ArgumentParser(description="Offset some set of columns by a fai.")
ap.add_argument("table", help="Table file")
ap.add_argument("fai", help="fai file")
ap.add_argument("--chromi", help="Column number of chromosome", type=int, default=0)
ap.add_argument("--posi", help="Column of position in chromosome", type=int, default=1)
args = ap.parse_args()

table = open(args.table)
fai = Tools.ReadFAIFile(args.fai)

for line in table:
    vals = line.split()
    chrom = vals[args.chromi]
    pos = int(vals[args.posi])
    offset = pos + fai[chrom][1]
    sys.stdout.write("\t".join(vals) + "\t" + str(offset) + "\n")
    
