#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import numpy as np
import pdb

ap = argparse.ArgumentParser(description="Look for replication timing transitions in the data.")
ap.add_argument("--input", help="bedGraph files with min 1kpb intervals", nargs="+", required=True)
ap.add_argument("--fai", help="input fai file for precomputed genome size", required=True)
ap.add_argument("--output", help="Output file.", default="/dev/stdout")
ap.add_argument("--window", help="Window to consider for transition.",default=10000)
ap.add_argument("--phase", help="Write phase bed to this file.", default=None)
ap.add_argument("--key", help="Order of phase input.", default=None)



args = ap.parse_args()

# First read in the fai file

fai = pd.read_csv(args.fai, sep='\t', header=None)
sys.stderr.write("Done reading fai.\n")

outputFile = open(args.output, 'w')

if (args.phase is not None):
    phaseFile = open(args.phase, 'w')

bedGraphFiles = args.input

if (args.key is not None):
    keys = args.key.split(',')
else:
    keys = [str(i) for i in range(0,len(bedGraphFiles))]
    
srcTables = [ pd.read_csv(filename, sep='\t', header=None) for filename in bedGraphFiles]
sys.stderr.write("Read source tables.\n")

tablesByChrom = [{} for t in srcTables]

for i in range(0,len(bedGraphFiles)):
    tablesByChrom[i] = { r[1][0] : np.zeros(r[1][1]/1000 + 1)  for r in fai.iterrows() }
sys.stderr.write("Organized by chromosome.\n")

for i in range(0,len(bedGraphFiles)):
    sys.stderr.write("Binning " + bedGraphFiles[i] + "\n")
    for chrom in fai[0]:
        for r in srcTables[i].loc[srcTables[i][0] == chrom].iterrows():
            tablesByChrom[i][chrom][r[1][1]/1000:r[1][2]/1000] = r[1][3]

sys.stderr.write("Finished binning.\n")

nPhases = len(tablesByChrom)

maxByChrom = { chrom : [ np.array([tablesByChrom[i][chrom][j] for i in range(0,nPhases)]).argmax() for j in range(0,len(tablesByChrom[0][chrom] )) ]  for chrom in fai[0] }



for chrom in fai[0]:
    l = len(maxByChrom[chrom])
    maxRunLength = np.zeros(l)
    maxRunLength[0] = 1
    
    vect = maxByChrom[chrom]

    for i in range(1,l):
        if (vect[i] == vect[i-1]):
            maxRunLength[i] = maxRunLength[i-1] + 1
        else:
            maxRunLength[i] = 1

    maxRun = maxRunLength[-1]
    
    for i in range(l,1,-1):
        cur = maxRunLength[i-1]
        maxRunLength[i-1] = maxRun
        if (cur == 1):
            maxRun = maxRunLength[i-2]
            
    maxRunLength[0] = maxRun

    for i in range(1,len(vect)):
        if (vect[i] != vect[i-1] and maxRunLength[i-1] > 1 and maxRunLength[i] > 1 ):
            try:
                outputFile.write(chrom + "\t" + str(i*1000) + "\t" + str((i+1)*1000) + "\t" + str(vect[i-1]) + "\t" + str(vect[i]) + "\t" + keys[vect[i-1]] + "\t" + keys[vect[i]] + "\n")
            except (IndexError):
                print str(i) + " is out of range for " + str(len(keys))

    if (args.phase is not None):
        regionStart = 0
        for i in range(1,len(vect)):
            if (vect[i] != vect[regionStart]):
                phaseFile.write(chrom + "\t" + str(regionStart*1000) + "\t" + str(i*1000) + "\t" + keys[vect[i]] + "\n")
                regionStart = i
        lasti = len(vect)-1
        lastv = vect[-1]
        phaseFile.write(chrom + "\t" + str(regionStart*1000) + "\t" + str(lasti*1000) + "\t" + keys[lastv] + "\n")
                

outputFile.close()
        


