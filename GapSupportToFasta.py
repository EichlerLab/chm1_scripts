\#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq, SeqRecord
import sys
import argparse
ap = argparse.ArgumentParser(description="Print fasta from support table")
ap.add_argument("table", help="Input table.")
ap.add_argument("fasta", help="Output fasta")
ap.add_argument("--type", help="insertion|deletion|both", default="both")
args = ap.parse_args()
    
tableName = args.table;
fastaName = args.fasta;



tableFile = open(tableName, 'r')
fastaFile = open(fastaName, 'w')

for line in tableFile:
    vals = line.split()
    support = int(vals[0])
    if (args.type == "both" or args.type == vals[1]):
        title = vals[2] + ":" + vals[3] + "-" + vals[4] + "/" + vals[1] + "/" + vals[0]
        for i in range(5,5+support):
            localTitle = title + '/' + str(i - 5)
            seq = SeqRecord.SeqRecord(Seq.Seq(vals[i]), id=localTitle, name="", description="")
            SeqIO.write(seq,fastaFile,"fasta")

fastaFile.close()
tableFile.close()
        
        
