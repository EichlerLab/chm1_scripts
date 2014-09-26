#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq

import sys
handle = open(sys.argv[1])

for seq in SeqIO.parse(handle, "fasta"):
    nc = seq.seq.count("C") + seq.seq.count("c")
    na = seq.seq.count("A") + seq.seq.count("a")
    ng = seq.seq.count("G") + seq.seq.count("g")
    nt = seq.seq.count("T") + seq.seq.count("t")
    l = len(seq.seq)
    print "{}\t{:2.2f}\t{:2.2f}\t{}".format(seq.id,float(nc+ng)/l, float(nt+na)/l, l)
