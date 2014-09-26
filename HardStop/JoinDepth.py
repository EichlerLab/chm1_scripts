#!/usr/bin/env python

import sys
import Tools
a = open(sys.argv[1])
b = open(sys.argv[2])

cov = {}
for line in a:
    v = line.split()
    (chr, pstart, pend) = Tools.ParseRegionStr(v[0])
    start = (pstart / 500)*500
    end = start + 500
    v[0] = "{}:{}-{}".format(chr, start, end)
    cov[v[0]] = float(v[1])

for line in b:
    v = line.split()
    c = "NA"
    if (len(v) != 3):
        continue
    
    if (v[0] in cov):
        c = str(cov[v[0]])
        v.append(c)
        sys.stdout.write("\t".join(v) + "\n")
    
