#!/usr/bin/env python

import argparse
import sys
import Tools
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="Take an insertion bed file, and print the surrounding sequence.")
ap.add_argument("bed", help="Input BED file.  This should be the breakpoints of the inversion")
ap.add_argument("genome", help="Input genome with a .fai file.")
ap.add_argument("contexts", help="Output file.")
ap.add_argument("--onlyflank", help="Only print flanking sequences.", action='store_true', default=False)
ap.add_argument("--window", help="Amount to store on sides.", type=int, default=1000)

args = ap.parse_args()

bedFile = open(args.bed)
contextFile = open(args.contexts, 'w')

fai = Tools.ReadFAIFile(args.genome+".fai")

genomeFile = open(args.genome)


for line in bedFile:
    vals  = line.split()
    chrom = vals[0]
    start = int(vals[1])

    
    prefixCoords = (max(0, start - args.window), start)
    suffixCoords = (start, min(fai[chrom][0], start+args.window))

    prefix = Tools.ExtractSeq((chrom, prefixCoords[0],prefixCoords[1]), genomeFile, fai)
    suffix = Tools.ExtractSeq((chrom, suffixCoords[0],suffixCoords[1]), genomeFile, fai)
    seqname="/".join(vals[0:3])
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(prefix+suffix), id=seqname, name="",description=""), contextFile, "fasta")
        

    

    
