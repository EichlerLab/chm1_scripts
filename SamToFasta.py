#!/usr/bin/env python
import pysam
import sys
import argparse
import Tools
import pdb
from Bio import SeqIO 
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="Select the read extending most into a gap.")
ap.add_argument("--sam", help="Input sam file, or stdin.", default="/dev/stdin")
ap.add_argument("--fasta", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

samFile = pysam.Samfile( args.sam, 'r')
outFile = open(args.fasta, 'w')
for aln in samFile.fetch():
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(aln.seq), id=aln.qname, name="", description=""), outFile, "fasta")
outFile.close()
