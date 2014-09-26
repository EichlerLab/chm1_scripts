#!/usr/bin/env python

import pysam
import sys
import argparse
import Tools
import pdb
import random

from Bio import SeqIO 
from Bio import SeqRecord
from Bio import Seq
import os

ap = argparse.ArgumentParser(description="Select the read extending most into a gap.")
ap.add_argument("bam", help="Indexed  file of alignments.")
ap.add_argument("region", help="Region to align", nargs="+")
ap.add_argument("--out", help="Write to this file", default="stdout")
ap.add_argument("--minqv", help="Minimum mapping quliaty", default=30, type=int)
ap.add_argument("--primary", help="Only store primary values.", default=False, action='store_true')

ap.add_argument("--max", help="Maximum number of reads to extract", default=1000, type=int)
ap.add_argument("--subsample", help="Subsample by this portion", default=None, action='store_true')
ap.add_argument("--fastq", help="Write fastq format.", default=False, action='store_true')
ap.add_argument("--filter", help="Only print a read if it is in this list.", default=None)
ap.add_argument("--slop", help="Add this to the region.", default=0, type=int)

args = ap.parse_args()
outFile = sys.stdout
if (args.out != "stdout"):
    outFile = open(args.out, 'w')

if (os.path.exists(args.region[0])):
    regionFile = open(args.region[0])
    args.region = [line.strip() for line in regionFile.readlines()]

if (args.filter is not None):
    with open(args.filter) as filterFile:
        filt = { l.strip(): True for l in filterFile }


printed = {}


for regionStr in args.region:
    region = Tools.FormatRegion(regionStr)

    
    
    if (region is None):
        print "malformatted region " + ' '.join(regionStr)
        sys.exit(1)

    region = (region[0], region[1]-args.slop, region[2]+args.slop)
    
    samFile = pysam.Samfile(args.bam)
    nAlns = samFile.count(region[0], region[1], region[2])
    if (nAlns <= args.max or args.subsample):
        if (nAlns > args.max and args.subsample ):
            lengths = []
            index = 0
            tmpPrinted = {}
            for aln in samFile.fetch(region[0], region[1], region[2]):
                if (aln.mapq < args.minqv):
                    continue
                if (args.primary and aln.flag & 256 != 0):
                    continue
                if (aln.qname in tmpPrinted):
                    continue
                
                tmpPrinted[aln.qname] = True 
                lengths.append((aln.tlen, index))
                index +=1
            lengths.sort()
            lengths.reverse()
            lastIndex = min(len(lengths), args.max)
            sys.stderr.write("last index: " + str(lastIndex) + "\n")
            lengths = lengths[0:lastIndex]
            indices = {i[1] : True for i in lengths}
            
        index = 0

        for aln in samFile.fetch(region[0], region[1], region[2]):

            if (aln.mapq < args.minqv):
                continue
            if (args.primary and aln.flag & 256 != 0):
                continue
            if (args.filter is not None and aln.qname not in filt):
                    continue

            if (aln.qname not in printed):
                rname = samFile.getrname(aln.tid)
                readname = aln.qname + " " + rname + ":" + str(aln.pos) + "-" + str(aln.aend)
                if ( nAlns <= args.max or (nAlns > args.max and args.subsample and index in indices)):
                    if (args.fastq == False):
                        SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(aln.seq), id=readname, name="", description=""), outFile, "fasta")
                    else:
                        outFile.write("@" +  readname + "\n")
                        outFile.write(aln.seq + "\n")
                        outFile.write("+\n")
                        outFile.write(aln.qual + "\n")
                printed[aln.qname] = True
                index += 1
    else:
        sys.stderr.write( regionStr + " contains too many reads: " + str(nAlns) + "\n")
        sys.exit(1)

if (args.out != "stdout"):
    outFile.close()
