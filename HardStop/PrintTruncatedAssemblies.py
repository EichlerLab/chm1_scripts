#!/usr/bin/env python

import Tools
import sys
import argparse
import pdb

ap = argparse.ArgumentParser(description="Print alignments that are not full length.")
ap.add_argument("sam", help="Input sam file.")
ap.add_argument("--names", help="File with names of gap-free alignments.", default=None)
ap.add_argument("--regions", help="Check assemblies from these regions.", default=None)
ap.add_argument("-v", help="Print complement", action='store_true', default=False)
ap.add_argument("--fraction", help="Min fraction of contig length that must be aligned.", type=float, default=0.9)

ap.add_argument("--onTarget", help="Print alignments that overlap the target given in the query name.", default=False, action='store_true')

args = ap.parse_args()

samFile = open(args.sam)
names = None
if (args.names is not None):
    nameFile = open(args.names)
    names = {}
    for line in nameFile:
        vals = line.split()
        names[vals[0].strip()] = line.strip()

regions = None        
if (args.regions is not None):
    regionFile = open(args.regions)
    regions = {}
    for line in regionFile:
        line = line.strip()
        regions[line] = True


def Lookup(title, names, regions):
    if (names is None and regions is None):
        return True
    
    if (names is not None):
        if (title in names):
            return True
        else:
            return False
    if (regions is not None):
        valRegion = title.split("_")[0]
        if (valRegion in regions):
            return True
        else:
            return False



def OnTarget(region, samValues):
    start = int(samValues[3])
    end   = start  + int(samValues[8])
    chrom = samValues[2]
    if (chrom != region[0]):
        return False

    if (start <= region[1] and region[2] <= end):
        return True
    return False

for line in samFile:
    if (line[0] == '@'):
        continue
    vals = line.split()
    start = int(Tools.GetKV("XS:i:", vals))
    end   = int(Tools.GetKV("XE:i:", vals))
    (contigStart, contigEnd) = vals[0].split('/')[-1].split("_")
    contigStart = int(contigStart)
    contigEnd   = int(contigEnd)

    l = contigEnd - contigStart

    printed = False
    region = Tools.ParseRegionStr(vals[0].split("_")[0])

    if (Lookup(vals[0], names, regions) == False):
        continue

    if (args.onTarget == True and OnTarget(region, vals) == False):
        continue

    if ( float(end - start)/l < args.fraction):
        printed = True
        if (args.v == False):
            print vals[0] + "\t" + "{:2.2f}".format(float(end - start)/l)  
        
    if (args.v and printed == False):
        
            print vals[0] +  "\t" + "{:2.2f}".format(float(end - start)/l)  
