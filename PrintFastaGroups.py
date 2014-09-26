#!/usr/bin/env python

import argparse
import sys
import os
import subprocess

ap = argparse.ArgumentParser(description="Print groups of fasta files that add to a target size.")
ap.add_argument("fasta", help="Input fasta files.", nargs="+")
ap.add_argument("--base", help="Basename for output.", default="ref")
ap.add_argument("--size", help="Target size", type=int, default=2000000000)

args = ap.parse_args()
total = 0
curFile = 0
files = []
for fileName in args.fasta:
    total += os.path.getsize(fileName)
    files.append(fileName)
    if (total >= args.size):
        print "cat " + " ".join(files) + " > " + args.base + "." + str(curFile) + ".fasta"
        total = 0
        curFile +=1
        files = []
