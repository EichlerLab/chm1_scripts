#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="Translate a region/bed line to a region.")
ap.add_argument("infile", nargs="?", help="input or none for stdin", default=None)
ap.add_argument("--slop", help="Add/subtract this amount from each value.", default=0, type=int)
args = ap.parse_args()

if (args.infile is not None):
    fh = open(args.infile)
else:
    fh = sys.stdin

for line in fh:
    vals = line.split()
    if (vals[0].find(":") >= 0):
        chrom = vals[0].split(":")[0]
        [starts,ends] = vals[0].split(":")[1].split("-")[0:2]
        start = int(starts)
        end   = int(ends)
    else:
        chrom = vals[0]
        start = int(vals[1])
        end   = int(vals[2])
    start = max(0,start - args.slop)
    end  =  end + args.slop

    print(chrom + ":" + str(start) + "-" + str(end))
