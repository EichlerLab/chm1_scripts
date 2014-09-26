#!/usr/bin/env python

import sys
import os
import argparse
ap = argparse.ArgumentParser(description="Print a file list using a bunch of read ids barcode/read/pos, by cutting first field of sam output.")
ap.add_argument("titles", help="barcode/read/pos")
ap.add_argument("output", help="file_name zmw")
ap.add_argument("--dir", help="Use this directory for input files.", default=".")
ap.add_argument("--fofn", help="Write fofn.", default=None)
args = ap.parse_args()
inFile = open(args.titles)
outFile = open(args.output, 'w')

barcodes = {}

fofnFile = None
if (args.fofn is not None):
    fofnFile = open(args.fofn, 'w')

for line in inFile:
    vals = line.split()
    title = vals[0].split('/')
    name  = title[0]
    hole  = int(title[1])
    if name not in barcodes:
        barcodes[name] = {}
    barcodes[name][hole] = True

fileNames = {}
for name in barcodes.keys():
    zmws = barcodes[name].keys()
    zmws.sort()
    suffix = ""
    for zmw in zmws:
        if (zmw >= 0 and zmw <= 54493):
            suffix = ".1.bax.h5"
        elif (zmw >= 54494 and zmw <= 108987):
            suffix = ".2.bax.h5"
        else:
            suffix = ".3.bax.h5"
        fileName = args.dir + "/" + name + suffix
        if (os.path.exists(fileName) == False):
            suffix = ".bas.h5"
            fileName = args.dir + "/" + barcodes[name] + suffix
            if (os.path.exists(fileName) == False):
                print "ERROR! Could not find a file for " + name
                print "Check your file structure and make sure the input is linked to the current directory."
                sys.exit(1)
        fileNames[fileName] = True
        outFile.write(fileName + " " + str(zmw) + "\n")
outFile.close()

if (fofnFile is not None):
    for fileName in fileNames.keys():
        fofnFile.write(fileName + "\n")
    fofnFile.close()


    
