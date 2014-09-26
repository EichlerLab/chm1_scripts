#!/usr/bin/env python

from Bio import SeqIO 
from Bio import SeqRecord
from Bio import Seq
import subprocess
import Tools
import sys
import time
import tempfile
import os
import argparse
ap = argparse.ArgumentParser(description = "Wrapper to plot a query against a region.")
ap.add_argument("query", help="query file, only dotplot first read.")
ap.add_argument("target", help="Target file, must have faidx.")
ap.add_argument("region", help="Region in chr:start-end format")
ap.add_argument("--word", help="Word size", type=int, default=13)
ap.add_argument("--slop", help="Adjust target by this much.", type=int, default=0)
args = ap.parse_args()


def Plot(queryName, targetName, region, output):
    with tempfile.NamedTemporaryFile(dir=".") as tempTarget:
        command = "samtools faidx {} {} ".format(targetName, region)
    

        subprocess.call(command.split(), stdout=tempTarget)
        command = "dottup -bsequence {} -asequence {} -wordsize {} -graph png".format(queryName, tempTarget.name, args.word)
        subprocess.call(command.split())
        command = "mv dottup.1.png {}.png".format(output)
        subprocess.call(command.split())


r = Tools.ParseRegionStr(args.region)
regionTup = [r[0], r[1], r[2]]
regionTup[1] = max(0, regionTup[1] - args.slop)
regionTup[2] = regionTup[2] + args.slop

region = "{}:{}-{}".format(regionTup[0], regionTup[1], regionTup[2])

index = 1
queryFile= open(args.query)
for record in SeqIO.parse(queryFile, "fasta"):
    with tempfile.NamedTemporaryFile(dir=".", delete=False) as tempQuery:
        SeqIO.write(record, tempQuery, "fasta")
        tempQuery.close()
        Plot(tempQuery.name, args.target, region, str(index))
        os.remove(tempQuery.name)

    with tempfile.NamedTemporaryFile(dir=".", delete=False) as tempQuery:
        seqrc = record.seq.reverse_complement()
        recordrc = SeqRecord.SeqRecord(seq=seqrc, id=record.id, description="", name="")
        SeqIO.write(recordrc, tempQuery, "fasta")
        tempQuery.close()
        Plot(tempQuery.name, args.target, region, str(index) + ".rc")
        os.remove(tempQuery.name)        
    index += 1
