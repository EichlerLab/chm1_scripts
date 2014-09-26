#!/usr/bin/env python

from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

import sys
import Tools
import argparse
import os
import re

ap = argparse.ArgumentParser(description="Fetch fasta sequences from indexed files.")

ap.add_argument("reference", help="FAI indexed reference file.")
ap.add_argument("region", help="UCSC style region (commas ok), or BED", nargs="+")
ap.add_argument("--slop", help="Add this to each side of the region.", default=0, type=int)
ap.add_argument("--prefix", help="Add this to the beginning of each region.", default=0, type=int)
ap.add_argument("--suffix", help="Add this to the end of each region.", default=0, type=int)
ap.add_argument("--out", help="Write here, default is stdout.", default=None)


args = ap.parse_args()


region = Tools.FormatRegion(args.region)

regionStr = "{}:{}-{}".format(region[0], region[1], region[2])

faiFile = args.reference + ".fai"
if (os.path.exists(faiFile) == False):
    print "The fasta file must be referenced and have a .fai index."
    
fai = Tools.ReadFAIFile(faiFile)
reffd = file(args.reference, 'r')
args.prefix -= args.slop
args.suffix += args.slop

region = Tools.AddPreSuf(region, fai, args.prefix, args.suffix)
seq = Tools.ExtractSeq(region, reffd, fai)

if (args.out is None):
    outFile = sys.stdout
else:
    outFile = open(args.out, 'w')

regionStr ="{}:{}-{}".format(region[0], region[1], region[2])
rec = SeqRecord.SeqRecord(Seq.Seq(seq), id=regionStr, name="", description="")
SeqIO.write(rec, outFile, "fasta")
if (outFile != sys.stdout):
    outFile.close()
