#!/usr/bin/env python

import argparse
ap = argparse.ArgumentParser(description="count masked repeats.")
ap.add_argument("table", help="Repeat table.")
ap.add_argument("--fraction", help="Calls must be at least this fraction to be made.", type=float, default=0.51)
args = ap.parse_args()


table = open(args.table, 'r')

table.readline()
table.readline()
table.readline()

indels = dict([("insertion", {}), ("deletion", {})])

for line in table:
    vals = line.split()
    title = vals[4]
    titleVals = title.split('/')
    region = titleVals[0]
    inordel = titleVals[1]
    if (region not in indels[inordel]):
        indels[inordel][region] = {}
    if (vals[9] not in indels[inordel][region]):
        indels[inordel][region][vals[9]] = 0
    indels[inordel][region][vals[9]] += 1


for k in ("insertion", "deletion"):
    for region in indels[k].keys():
        total = 0.0;
        for op in indels[k][region].keys():
            total += indels[k][region][op]
        for op in indels[k][region].keys():
            frac = indels[k][region][op] / total
            if ( frac >= args.fraction):
                print k + "\t" + region + "\t" + op + "\t" + str(frac)
        
        

