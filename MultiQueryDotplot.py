#!/usr/bin/env python

from Bio import SeqIO
from Bio import SeqRecord
import subprocess
import tempfile
import argparse
import sys
ap = argparse.ArgumentParser(description="Run dotplots for a few sequences.")

ap.add_argument("query", help="Multi-fasta query file.")
ap.add_argument("target", help="target.")

args = ap.parse_args()

tmpQueryName = tempfile.mktemp(suffix=".query.fasta", dir=".")

h = open(args.query)
index = 0
for seq in SeqIO.parse(h, "fasta"):
    tmpFile = open(tmpQueryName, 'w')
    SeqIO.write(seq, tmpFile, "fasta")
    tmpFile.close()

    command = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DotPlot.py --query {} --target {} --matches dot:11 --savefig dotplot.{}.png".format(tmpQueryName, args.target, index)
    print command
    subprocess.call(command.split())
    index += 1


