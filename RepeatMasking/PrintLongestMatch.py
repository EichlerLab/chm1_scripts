#!/usr/bin/env python

import argparse
ap = argparse.ArgumentParser(description="Count the repeats that are aligned.")
ap.add_argument("input", help="file of file names (.m1 alignments)")
args = ap.parse_args()

inFile = open(args.input)

for fileName in inFile:
    m1File= open(fileName.strip())
    alns = {}    
    for line in m1File:
        vals = line.split()
        name = "/".join(vals[0].split("/")[0:2])
        alnLength = int(vals[10]) - int(vals[9])
        totalLength = int(vals[11])
        target = vals[1]
        ident  = vals[5]
        if (name not in alns):
            alns[name] = [target, ident, alnLength, totalLength, float(alnLength)/totalLength]
        else:
            if (alns[name][2] < alnLength):
                alns[name] = [target, ident, alnLength, totalLength, float(alnLength)/totalLength]

        for aln,vals in alns.iteritems():
            print aln + "\t" + "\t".join([str(i) for i in vals])
                
        
