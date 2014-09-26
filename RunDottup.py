#!/usr/bin/env python

from Bio import SeqIO
import sys
import subprocess
import tempfile
import argparse
ap = argparse.ArgumentParser(description="Self dot of everything in a fasta file.")
ap.add_argument("fasta", help="Input multifasta file.")
ap.add_argument("--base", help="Output base", default="dottup.")
ap.add_argument("--target", help="Use this target instead of self.", default=None)
ap.add_argument("--n", help="Only at most this many", default=0, type=int)
args = ap.parse_args()
fastaFile = open(args.fasta)

index = 0
for rec in SeqIO.parse(fastaFile, "fasta"):
    tf = tempfile.NamedTemporaryFile(delete=False)
    SeqIO.write(rec, tf, "fasta")
    tf.close()
    target = tf.name

    if (args.target is not None):
        target = args.target
    print tf.name
    print target
    print rec.id
    tmpGraph = tempfile.NamedTemporaryFile();
    command = "dottup -asequence {} -bsequence {} -wordsize 11 -graph ps -goutfile {}".format(tf.name, target, tmpGraph.name)
    subprocess.call(command.split())
    outName = rec.id
    outName = outName.replace(":", "_")
    outName = outName.replace("/", "_")

    command = "convert -density 300 -rotate 90 -flatten -background white {}.ps {}".format(tmpGraph.name, args.base +outName + ".png")
    subprocess.call(command.split())

    command = "rm -f {}".format(tf.name)
    subprocess.call(command.split())    
    index += 1
    if (index == args.n):
        break
