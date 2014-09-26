#!/usr/bin/env python

import sys
import Tools
import argparse
ap = argparse.ArgumentParser(description="Print tabular version of alignments with intervals.")
ap.add_argument("sam", help="Input sam files", nargs="+")
args = ap.parse_args()

vals = []
for samFile in args.sam:
    sys.stderr.write(samFile + "\n")
    alignFile = open(samFile)

    for line in alignFile:
        if (line[0] == "@"):
            continue
        sam = Tools.SAMEntry(line)
        vals.append((sam.title, sam.flag, sam.qStart, sam.qEnd, sam.tName, sam.tStart, sam.tEnd))
#        print "\t".join([str(i) for i in (sam.title, sam.flag, sam.qStart, sam.qEnd, sam.tName, sam.tStart, sam.tEnd)])

vals.sort

for v in vals:
    print "\t".join(str(i) for i in v)

