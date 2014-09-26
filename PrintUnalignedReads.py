#!/usr/bin/env python

import sys
from Bio import SeqIO
import os
import argparse
ap = argparse.ArgumentParser(description="Print reads with titles not matching the unaligned file.")
ap.add_argument("key", help="The first field is the read name. File of read names.")
ap.add_argument("fasta", help="Fasta input", nargs="+")
ap.add_argument("--output", help="Output file, default is input.unmapped.")

args= ap.parse_args()

namesFile = open(args.key)
names = {}
i=0
for line in namesFile:
    vals = line.split()
    names[vals[0]] = True
    i+=1
    if (i % 100000 == 0):
        print "Indexed " + str(i)

if (args.output is not None):
    unmappedReadsFile = open(args.output, 'w')
    
for readsFileName in args.fasta:
    readsFile = open(readsFileName)
  
    localFile = os.path.basename(readsFileName)
    if (args.output is None):
        unmappedReadsFile = open(localFile  + ".unmapped", "w")
    for read in SeqIO.parse(readsFile, "fasta"):
        if (read.id not in names):
            SeqIO.write(read, unmappedReadsFile, "fasta")
    if (args.output is None):
        unmappedReadsFile.close()



    sys.stderr.write(readsFileName + "\n")
    
