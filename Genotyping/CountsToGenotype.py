#!/usr/bin/env python


import sys
import numpy as np
import argparse

ap = argparse.ArgumentParser(description="print mean or median of counts.")
ap.add_argument("counts", help="Input file in format 'label count', where label happens many times")
ap.add_argument("out", help="Print genotype to this file.")
ap.add_argument("--stat", help="Statistic to consider. (mean, median)", default="mean")
ap.add_argument("--max", help="Exclude entries with a value greater than this max.", default=None, type=int)
ap.add_argument("--coverage", help="Coverage table for source data.", default=1.0, type=float)
ap.add_argument("--minValue", help="Count queries with at least this many reads towards average.", default=0, type=int)
args = ap.parse_args()

countsFile = open(args.counts)
outFile = open(args.out, 'w')
labels = []
counts = []
for line in countsFile:
    vals = line.split()
    if (args.max is None or int(vals[2]) < args.max):
        labels.append(vals[0])
        counts.append(int(vals[2]))

indices = [0]
counts = np.array(counts)

for i in range(len(labels)-1):
    if (labels[i] != labels[i+1]):
        indices.append(i+1)
        

indices.append(len(labels))


def GetMean(values, minVal):
    vals = np.array(values[np.where(values >= minVal)])/args.coverage
    if (np.size(vals) == 0):
        return 0
    else:
        return np.mean(vals)

def GetMedian(values, minVal):
    vals = np.array(values[np.where(values >= minVal)])/args.coverage
    if (np.size(vals) == 0):
        return 0
    else:
        return np.median(vals)

    
if (args.stat == "mean"):
    
    stats = [ GetMean(counts[indices[i]:indices[i+1]], args.minValue) for i in range(0,len(indices)-1)]
elif (args.stat == "median"):
    stats = [ GetMedian(counts[indices[i]:indices[i+1]], args.minValue) for i in range(0,len(indices)-1)]
for i in range(0,len(indices)-1):
    outFile.write("{} {:d} {:d} {}\n".format(labels[indices[i]], int(stats[i]), int(indices[i+1]-indices[i]), i+1))
outFile.close()
    



    
