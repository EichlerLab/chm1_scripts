#!/usr/bin/env python

import Tools
import sys

import argparse
import numpy as np


ap = argparse.ArgumentParser(description="Given a table of query sequence names and coordinates in the query sequence, test for presence of that query in alignments to another species.")
ap.add_argument("queryTable", help="Table of query sequence names, start and end coordinate (bed format)")
ap.add_argument("alignedQueries", help="SAM file of alignment queries, not necessarily in the same order as table.")
ap.add_argument("output", help="Output file name.")
ap.add_argument("--bedplus", help="Input is a bed file, with the value of this as the value of the bed.", type=int, default=None)

args = ap.parse_args()


queryTable = open(args.queryTable)
samFile = open(args.alignedQueries)
outputFile = open(args.output, 'w')


queries = {}

if (args.bedplus is None):
    for line in queryTable:
        vals = line.strip().split()
        queries[vals[0]] = vals
else:
    #
    # The query table is instead a BED file that was used to generate
    # the alignments.  The value of bedplus should be the same as the
    # value used to generate the fasta file of sequence plus flanks.
    #
    for line in queryTable:
        vals = line.split()
        queryName = "/".join(vals[0:3])
        res = [queryName, args.bedplus, args.bedplus + int(vals[2]) - int(vals[1])]
        queries[queryName] = res
        
keys = queries.keys()

for alnLine in samFile:
    if (alnLine[0] == '@'):
        continue
    aln = Tools.SAMEntry(alnLine)
    alnTitle = "/".join(aln.title.split("/")[0:3])
    if (alnTitle in queries):
        queryAln = np.zeros(len(aln.seq))
        qPos = aln.qStart

        for i in range(0,len(aln.ops)):
            if (aln.ops[i] == 'S' or aln.ops[i] == 'H'):
                continue
            if (aln.ops[i] == 'M'):
                queryAln[qPos:qPos + aln.lengths[i]] = 1
                qPos += aln.lengths[i]
            if (aln.ops[i] == 'I'):
                qPos += aln.lengths[i]

        queryFocusStart = int(queries[alnTitle][1])
        queryFocusEnd = int(queries[alnTitle][2])
        
        coverageSum = sum(queryAln[queryFocusStart:queryFocusEnd])

        queries[alnTitle].append( coverageSum / (queryFocusEnd - queryFocusStart))

        
for k in keys:
    if (len(queries[k]) == 3):
        queries[k].append(0)
    outputFile.write("\t".join([str(i) for i in queries[k]]) + "\n")
