#!/usr/bin/env python

import os
import sys
import argparse
import Tools

ap = argparse.ArgumentParser(description="Find an alignment file.")
ap.add_argument("filelist", help="Barcode of read, long string staring with m1....", nargs="?", default=None)
ap.add_argument("--dir", help="Look here for file.", default=".")
ap.add_argument("--findex", help="Write in findex format for input into writeHDFSubset." ,default=False,action='store_true')
ap.add_argument("--dotplot", help="Write in dotplot parameter format.", default=False, action='store_true')
args = ap.parse_args()


fileList = []

if (args.filelist is None):
    for line in sys.stdin.readlines():
        fileList.append(line)

else:    
    for fl in args.filelist:
        if (os.path.exists(fl)):
            print "this exists: " + fl
            fileListFile = open(fl, 'r')
            fileList.extend([line.strip() for line in fileListFile.readlines()])
        else:
            if (fl.find(';') >= 0):
                fileList.extend(fl.split(';'))
            else:
                fileList.append(fl)

existingReads = {}
dotplotArgs = []
for line in fileList:
    (barcode,zmw,coords) = Tools.ParseReadTitle(line.split()[0])
    read = barcode+ "/" + str(zmw)
    if read in existingReads.keys():
        continue
    else:
        existingReads[read] = True

    filename = Tools.FindBasFile(barcode, zmw, args.dir)

    if (filename is None):
        print "Did not find " + read
        sys.exit(1)
    else:
        if (args.findex == True):
            print filename + " " + str(zmw)
        elif (args.dotplot == True):
            dotplotArgs.append("bas:{}:{}".format(filename, zmw))
        else:
            print filename

if (args.dotplot == True):
    print " ".join(dotplotArgs)

    

