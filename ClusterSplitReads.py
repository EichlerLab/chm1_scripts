#!/usr/bin/env python

import bisect
import argparse
import sys
import pdb

def FindInterval(chrom, start, clusters):
    s = (start, start, 0)
    if (chrom not in clusters):
        return None
    sb = bisect.bisect_left(clusters[chrom], s)
#    pdb.set_trace()
    while (sb < len(clusters[chrom]) and sb > 0 and clusters[chrom][sb][0] > start):
        sb -= 1
    if (sb >= len(clusters[chrom])):
        return None
    if (clusters[chrom][sb][0] <= start and clusters[chrom][sb][1] >= start):
        return clusters[chrom][sb][2]
    return None



def ReadIntervals(clusterFileName):
    fh = open(clusterFileName)
    intervals = {}
    i = 0
    allIntervals = []
    for line in fh:
        vals = line.split()
        if (vals[0] not in intervals):
            intervals[vals[0]] = []
        intervals[vals[0]].append((int(vals[1]), int(vals[2]), i))
        allIntervals.append((vals[0], vals[1], vals[2]))
        i += 1
    return (intervals, allIntervals)


ap = argparse.ArgumentParser(description="Cluster breakpoints hit by split reads.")
ap.add_argument("table", help="bed file of split read coordinates.  Reads should be paired, ll, or rr, or rl. Fields 1-3 are bed coordinates, field 4 is the read name, field 5 is length of hard-stop, and 6 is left/right.")
ap.add_argument("intervals", help="intervals")

ap.add_argument("--min", help="Only print clusters with at least this count.", default=5,type=int)
args = ap.parse_args()

(intervalHash, intervals)   = ReadIntervals(args.intervals)

table = open(args.table)
clusters = {}
lines = table.readlines()
i = 0
while (i < len(lines)):
    vals1 = lines[i].split()
    vals2 = lines[i+1].split()
    if (vals1[3] != vals2[3]):
        i+=1
        print vals1[3] + " is not paired"
        continue
    

    i1 = FindInterval(vals1[0], int((int(vals1[1]) + int(vals1[2]))/2), intervalHash)
    i2 = FindInterval(vals2[0], int((int(vals2[1]) + int(vals2[2]))/2), intervalHash)
#    print str(i1) + " " + str(i2)
    if (i1 is not None and i2 is not None):
        if (i1 not in clusters):
            clusters[i1] = {}
        if (i2 not in clusters[i1]):
            clusters[i1][i2] = []
        clusters[i1][i2].append(vals1[3])
    i += 2
    
for i1 in clusters.keys():
    for i2 in clusters[i1].keys():
        if (i1 == i2):
            continue
        if (len(clusters[i1][i2]) > args.min):
            print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(intervals[i1][0], intervals[i1][1], intervals[i1][2], intervals[i2][0], intervals[i2][1], intervals[i2][2], len(clusters[i1][i2]))
