#!/usr/bin/env python

import sys

if (len(sys.argv) != 4):
    print "usage: AnnotateGapBedWithCensorMap.py gaps.bed gaps.file gaps.out.bed "
    print "fields 1-6 are: chrom<\t>start<\t>end<\t>insertion/deletion<\t>length<\t>seq."
    print "the output file has the 7th field as the annotation, 8th field as the remainder, and 9th field as the percent repeat."
    
gapBed = open(sys.argv[1])
mapFile = open(sys.argv[2])
gapBedOut = open(sys.argv[3], 'w')

entries = {}
order = []
for line in gapBed:
    vals = line.split()
    name = "/".join(vals[0:3])
    entries[name] = vals
    entries[name].append([])
    order.append(name)

for line in mapFile:
    vals = line.split()

    if (vals[0] not in entries):
        print "name: " + vals[0] + " not found "
        continue
    name = vals[0]
    last = len(entries[name]) - 1

    entries[name][last].append(vals[3])
    repStart = int(vals[1])-1
    repEnd   = int(vals[2])
    substr = entries[name][5][repStart:repEnd]
    substr = substr.lower()
    prefix= entries[name][5][0:repStart]
    suffix=entries[name][5][repEnd:]
    entries[name][5] = prefix+substr+suffix


for n in order:
    if (len(entries[n][len(entries[n])-1]) == 0):
        entries[n][len(entries[n])-1].append("NONE")

    seq = entries[n][5]
    nLower = seq.count("a") + seq.count("c") + seq.count("t") + seq.count("g")
    frac = float(nLower) / len(seq)
    gapBedOut.write("\t".join(entries[n][0:6]) +"\t" + ";".join(entries[n][len(entries[n])-1]) + "\t" + "\t".join(entries[n][6:8]) + "\t" + "{:2.2f}".format(frac) + "\n")
