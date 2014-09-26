#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="Annotate a gap bed file with an associated TRF table.")
ap.add_argument("bed", help="Input BED file.")
ap.add_argument("trf", help="Input TRF table file.")
ap.add_argument("bed_out", help="Output BED file.")
args = ap.parse_args()

trfFile = open(args.trf, 'r')
prevSeqName = ""

annotations = {}
for line in trfFile:
    if (line[0] == '@'):
        seqName = line[1:].strip()
    else:
        vals = line.split()
        if (seqName not in annotations):
            annotations[seqName] = (int(vals[1]) - int(vals[0]), vals[13])
        else:
            gapLen = int(vals[1]) - int(vals[0])
            if (annotations[seqName][0] < gapLen):
                annotations[seqName] = (gapLen, vals[13])
print "done processing annotations"

bedFile = open(args.bed, 'r')
bedOutFile = open(args.bed_out, 'w')

for line in bedFile:
    vals = line.split()
    seqTitle = '/'.join(vals[0:3])
    if (seqTitle in annotations):
        bedOutFile.write('\t'.join(vals) + "\t" + str(annotations[seqTitle][0]) + "\t" + annotations[seqTitle][1] + "\n")
    else:
        bedOutFile.write('\t'.join(vals) + "\t0\tNO_TRF\n")
        
