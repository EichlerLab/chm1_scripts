#!/usr/bin/env python

import pysam
import sys
import argparse
import Tools

ap = argparse.ArgumentParser(description="Print coordinates of possibly only one subread.")
ap.add_argument("sam", help="Input sam file.")
ap.add_argument("bed", help="Output bed file.")
ap.add_argument("--best", help="Output only the best alignment", action='store_true', default=False)
ap.add_argument("--features", help="Additional features to annotate BED line", nargs="+", default=[])
ap.add_argument("--minqv", help="Minimum mapping quality", type=int, default=0)
args = ap.parse_args()

sf = open( args.sam, "r" )

if (args.bed != "stdout"):
    bedF = open(args.bed, 'w')
else:
    bedF = sys.stdout

prevBase = ""

def GetTLen(cigar):
    tLen = 0
    for op in cigar:
        if (op[0] == 0 or op[0] == 1):
            tLen += op[1]
    return tLen


prevBase = ""
featureStr = '\t'.join(args.features)

for line in sf:
    if (line[0] == '@'):
        continue
    aln = Tools.SAMEntry(line)
    name = aln.title
    base = '/'.join(name.split("/")[0:2])
    if (aln.mapqv < args.minqv):
        continue
    qStart = Tools.GetKV("XS:i:", aln.vals[11:])
    qEnd   = Tools.GetKV("XE:i:", aln.vals[11:])

    bedF.write(aln.tName + "\t" + str(aln.tPos)  + "\t" + str(int(aln.tPos) + int(aln.vals[8])) + "\t" + aln.title + "\t" + str(aln.strand) + "\t" + str(aln.qStart) + "\t" + str(aln.qEnd) + "\t" + str(int(qEnd) - int(qStart)) + "\n")


if (bedF != sys.stdout):
    bedF.close()
            

    
