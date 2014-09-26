#!/usr/bin/env python

import sys

import argparse
import pdb

def OnTarget(onTarget, vals):
    if (onTarget is None):
        return True
    start = int(vals[3])
    chrom = vals[2]
    end   = start + int(vals[8])
    for i in range(0,len(onTarget)):
        if (onTarget[i][0] == chrom and (onTarget[i][1] > start and onTarget[i][2] < end) or (onTarget[i][1] < start and onTarget[i][2] > start) or (onTarget[i][1] < end and onTarget[i][2] >= end)):
            return True
    return False
        
ap = argparse.ArgumentParser()
ap.add_argument("sam", help="alignments.")
ap.add_argument("--names", help="Print with these names.", default=None)
ap.add_argument("--regions", help="Filter region from name, and check against this list.", default=None)
ap.add_argument("-v", help="Opposite of logic.", action='store_true', default=False)
ap.add_argument("--onTarget", help="Print alignments if they are in this target list.", default=None)


args = ap.parse_args()
names = None
if (args.names is not None):
    name = open(args.names)
    names = {}
    for line in name:
        names[line.split()[0]] = True

regions = None
if (args.regions is not None):
    rf = open(args.regions)
    regions = {}
    for line in rf:
        vals = line.split()
        rn = vals[0]
        regions[rn] = True

onTarget = None
if (args.onTarget is not None):
    tf = open(args.onTarget)
    onTarget = []
    for line in tf:
        vals = line.split()
        onTarget.append((vals[0], int(vals[1]), int(vals[2])))
    
sam = open(args.sam)
for line in sam:
    if (line[0] == '@'):
        sys.stdout.write( line)
        continue
    
    doPrint = True
    vals = line.split()
    if (names is not None and vals[0] not in names):
        doPrint = False

    if (regions is not None):
        rgn = vals[0].split("_")[0]
        if (rgn not in regions):
            doPrint = False

    if (OnTarget(onTarget, vals) == False):
        doPrint = False
        
    if (doPrint and args.v == False):
        sys.stdout.write(line)

    if (doPrint == False and args.v == True):
        sys.stdout.write(line)
