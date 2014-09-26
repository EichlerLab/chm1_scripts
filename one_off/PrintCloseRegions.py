#!/usr/bin/env python


import sys
import argparse
ap = argparse.ArgumentParser(description="Join bed lines separated by no more than N kb.")

fh = open(sys.argv[1])
intvStart = None
for line in fh:
    vals = line.split()
    (chrom, start, end) = (vals[0], int(vals[1]), int(vals[2]))
    if (intvStart is not None and chrom == intvStart[0] and 
    
