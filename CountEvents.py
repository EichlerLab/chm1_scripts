#!/usr/bin/env python

import argparse
import sys
import Tools
import numpy as np

ap = argparse.ArgumentParser(description="Print the counts of events binned in the genome, one event per bed line.")
ap.add_argument("bed", help="Input bed file.")
ap.add_argument("genome", help="Genome index file.")
ap.add_argument("--window", help="Window to bin.", default=10000)



args = ap.parse_args()

fai = Tools.ReadFAIFile(args.genome)

def NBins(length, window):
    return length / window + length % window

def GetBin(chrom, chrommap, start, window):
    return (chrommap[chrom], start/window)
           
bins = [ np.zeros(NBins(fai[i][0], args.window)) for i in fai.keys() ]
chrommap = { fai.keys()[i]: i for i in range(len(fai.keys())) }
          


bedFile = open(args.bed)
for line in bedFile:
    vals = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end   = int(vals[2])
    (i, pos) = GetBin(chrom, chrommap, start, args.window )
    bins[i][pos] += 1

for chrom in fai.keys():
    chromi = chrommap[chrom]
    for i in range(len(bins[chromi])):
        if (bins[chromi][i] > 0):
            print chrom + "\t" + str(i*args.window) + "\t" + str(bins[chromi][i])
