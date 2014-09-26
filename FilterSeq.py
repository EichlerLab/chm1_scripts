#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import re

ap = argparse.ArgumentParser(description="Fetch fasta sequences from indexed files.")

ap.add_argument("input", help="input sequences")
ap.add_argument("--names", help="ids of reads", default = None)
ap.add_argument("--regions", help="ids contain regions", default=None)
ap.add_argument("--out", help="output file", default="/dev/stdout")

args = ap.parse_args()

inFile = open(args.input)
outFile = open(args.out, "w")
names = None
if (args.names is not None):
    names = {}
    a = open(args.names)
    for line in a:
        names[line.split()[0]] = True
    a.close()

regions =None
regionre = re.compile("(chr.*:\d+\-\d+).*")
if (args.regions is not None):
    regions = {}
    r = open(args.regions)
    for line in r:
        regions[line.split()[0]] = True
    r.close()
    

for record in SeqIO.parse(inFile, "fasta") :
    doPrint = False
    if (names is not None and record.id in names):
        doPrint = True

    if (regions is not None):
        m = regionre.match(record.id)
        if (m is not None and len(m.groups()) > 0):
            rgn = m.groups()[0]
            if (rgn in regions):
                doPrint = True
            
    if (doPrint):
        SeqIO.write(record, outFile, "fasta")
        
