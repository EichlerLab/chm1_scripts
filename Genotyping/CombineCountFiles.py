#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Combine multiple .count files into one.")
ap.add_argument("--input", help="input files", nargs="+")
ap.add_argument("--output", help="output files")

args = ap.parse_args()


def GetLines(filename):
    f = open(filename)
    return f.readlines()

inTables = [ GetLines(name) for name in args.input ]

outFile = open(args.output, 'w')

for i in range(0,len(inTables[1])):
    total = 0
    for j in range(len(inTables)):
        vals = inTables[j][i].split()
        total += int(vals[2])
    outFile.write(vals[0] + "\t" + vals[1] + "\t" + str(total) + "\t" + vals[3] + "\n")

outFile.close()
