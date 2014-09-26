#!/usr/bin/env python

from Bio import SeqIO
from Bio import Restriction
import argparse
import sys
import numpy

ap = argparse.ArgumentParser(description="Run an in-silico digest of a set of sequences.")
ap.add_argument("--motifs", help="Restriction site pattern (exact)", nargs="+")
ap.add_argument("--enzymes", help="Enzymes to search (eg. DraI, SwaI),", nargs="+")
ap.add_argument("--contigs", help="Input files (1 or more).", nargs="+")

opts = ap.parse_args()
inFileIndex = 0

def CountGC(seq, start, end):
    ng = seq.count("G", start, end) + seq.count("g", start, end)
    nc = seq.count("C", start, end) + seq.count("c", start, end)
    return ng + nc

for inFileName in opts.contigs:
    inFile = open(inFileName)
    contigSites = [ [] for m in opts.enzymes ]
    nSites = {}
    for contig in SeqIO.parse(inFile, "fasta"):
        sys.stderr.write("Done reading " + inFileName + "\n")
        mi = 0
        
        for e in opts.enzymes:
            enzyme = getattr(Restriction, e)
            contigSites[mi] = enzyme.search(contig.seq)
            mi += 1
        
        if (len(opts.enzymes) > 1):
            for i in range(1,len(opts.enzymes)):
                u = numpy.union1d(numpy.array(contigSites[0]),numpy.array(contigSites[i]))
            contigSites[0] = u
        for i in range(0,len(contigSites[0])-1):
            if (contigSites[0][i] in nSites or contigSites[0][i+1] in nSites):
                continue
            cutLen = contigSites[0][i+1] - contigSites[0][i]
            if (cutLen > 2):
                gc = CountGC(contig.seq, contigSites[0][i], contigSites[0][i+1])
                n  = contig.seq.count("N", contigSites[0][i], contigSites[0][i+1])
                if ((float(n)/cutLen) > 0.25):
                    continue
                print str(inFileIndex) + " " + str(contigSites[0][i]) + " " + str(cutLen) + " " + str(gc) + " " + str("{:2.2f}".format(float(gc)/cutLen)) + " " + str(n) + " " + str("{:2.2f}".format(float(n)/cutLen))
                

        inFileIndex += 1


