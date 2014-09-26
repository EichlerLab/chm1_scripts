#!/usr/bin/env python

from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import Tools
import sys
import os
import glob
import argparse
ap = argparse.ArgumentParser(description="Concatenate Celera assemblies into one output file.")
ap.add_argument("--base", help="Basename of table.  For each entry in the table, will search for a contig file.", default="rgn_")
ap.add_argument("--out", help="Output file name", default="stdout")
ap.add_argument("--fasta", help="Fasta file name", default="region/region.ctg.consensus.fasta")
ap.add_argument("--maxNumContigs", help="Only print assemblies if they have this number of contigs or fewer.", default=None, type=int)

args = ap.parse_args()
outFile = sys.stdout

if (args.out != "stdout"):
    outFile = open(args.out, 'w')
    


index = 1
regionSearch = args.base + "*"
regionDirs = glob.glob(regionSearch)
for regionDir in regionDirs:
    if (os.path.isdir(regionDir) == False):
        continue
    assemblyFile = regionDir + "/" + args.fasta
    regionFile = regionDir + "/region.txt"
    rf = open(regionFile)
    region = rf.readline()
    region = region.strip()
    index += 1
    if (os.path.exists(assemblyFile) == False):
        sys.stderr.write(assemblyFile + " not found.\n")
        continue

    asmFile = open(assemblyFile, 'r')
    seqs = [ seqrec for seqrec in SeqIO.parse(asmFile, "fasta") ]
    if (args.maxNumContigs is None or len(seqs) <= args.maxNumContigs):
        for i in range(0,len(seqs)):
            seqs[i].id = region + "_" + seqrec.id
            seqs[i].name = ""
            seqs[i].description = ""
            SeqIO.write(seqs[i], outFile, "fasta")
            
if (args.out != "stdout"):
    outFile.close()
