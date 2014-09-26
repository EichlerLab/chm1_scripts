#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

import Tools

ap = argparse.ArgumentParser(description="Print fasta from a sam file, possibly manipulating sequences.")
ap.add_argument("--sam", help="Input sam file", default="/dev/stdin")
ap.add_argument("--maskAligned", help="Mask the portion of the read that is aligned.", default=False, action='store_true')
ap.add_argument("--out", help="Output file", default="/dev/stdout")

args = ap.parse_args()

sam = open(args.sam)
fasta = open(args.out, 'w')
for line in sam:
    if (len(line) > 0 and line[0] == '@'):
        continue
    entry = Tools.SAMEntry(line)

    if (args.maskAligned):
        alignStart = 0
        alignEnd = len(entry.seq)
        i = 0
        while (i < len(entry.ops)):
            if (entry.ops[i] == 'S'):
                alignStart += entry.lengths[i]
                break
            i+=1
        j = len(entry.ops) - 1
        while (j > i):
            if (entry.ops[j] == 'S'):
                alignEnd -= entry.lengths[i]
                break
            j-=1

        newSeq = entry.seq[0:alignStart] + 'N'*(alignEnd - alignStart) + entry.seq[alignEnd:]
        entry.seq = newSeq
    rec = SeqRecord.SeqRecord(Seq.Seq(newSeq), id=entry.title, name="", description="")
    SeqIO.write(rec, fasta, "fasta")

fasta.close()
            
    
