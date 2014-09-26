#!/usr/bin/env python

import Tools
import sys
from Bio import SeqIO


seqIn = open(sys.argv[1])
seqOut = open(sys.argv[2], 'w')
regionStr = sys.argv[3]


region = Tools.ParseRegionStr(regionStr)
for record in SeqIO.parse(seqIn, "fasta"):
    
    record.seq = record.seq[0:region[1]]+ 'N'*(region[2]-region[1])  + record.seq[region[2]:]
    SeqIO.write(record, seqOut, "fasta")
    break





