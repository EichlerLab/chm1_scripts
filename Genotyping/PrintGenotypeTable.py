#!/usr/bin/env python

import argparse
import sys
import numpy as np

ap = argparse.ArgumentParser(description="Given a bunch of count files, print out by row, and maybe some descriptive stats.")
ap.add_argument("input", help="Input files.")
ap.add_argument("output", help="Output files.")
ap.add_argument("--statsBySample", help="Print statistics of nonzero coverage by sample.", action='store_true', default=False)
ap.add_argument("--statsByLocus", help="Print statics of loci.", action='store_true', default=False)


args = ap.parse_args()

output = open(args.output, 'w')

def FileToRow(fileName):
    f = open(fileName)
    l = f.readlines()
    v = [int(i.split()[1]) for i in l]
    return (fileName, v)

inputFile = open(args.input)
genotypeFiles = inputFile.readlines()


genotypes = [FileToRow(i.strip()) for i in genotypeFiles]

if (args.statsBySample == True):
    for i in range(0,len(genotypes)):
        total = 0
        n = 0
        nonzero = []
        for j in range(len(genotypes[i][1])):
            if (genotypes[i][1][j] != 0):
                nonzero.append(genotypes[i][1][j])
        a = np.asarray(nonzero)
        output.write("{} {:2.2f} {:2.2f}\n".format(genotypes[i][0], np.mean(a), np.std(a)))
    output.close()
    sys.exit(0)

if (args.statsByLocus == True):
    for i in range(0,len(genotypes[0][1])):
        total = 0
        n = 0
        nonzero = []
        for j in range(len(genotypes)):
            nonzero.append(genotypes[j][1][i])
        a = np.asarray(nonzero)
        output.write("{} {:2.2f} {:2.2f}\n".format(i, np.mean(a), np.std(a)))
    output.close()
    sys.exit(0)
        

for i in range(0,len(genotypes)):
    output.write(genotypes[i][0] + "\t" + "\t".join([str(v) for v in genotypes[i][1]]) + "\n")
output.close()
