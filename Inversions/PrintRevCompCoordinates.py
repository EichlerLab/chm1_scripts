#!/usr/bin/env python

import subprocess
import sys
import argparse

ap = argparse.ArgumentParser(description="Print inversion coordinates, assembly coordinates, or inversion lengths.")
ap.add_argument("alignfile", help=".m1 file of assemblies with inversions")
ap.add_argument("--length", help="Print length of inversion", default=False, action='store_true')
ap.add_argument("--asm", help="Print inversion coordinates in assembly.", default=False, action='store_true')
ap.add_argument("--gen", help="Print inversion coordiantes in genome.", default=False, action='store_true')
ap.add_argument("--asmgen", help="Print assembly coordinates in genome.", default=False, action='store_true')
ap.add_argument("--first", help="Print first half of inversion coordinates in assembly.", default=False, action='store_true')
ap.add_argument("--second", help="Print second half of inversion coordinates in assembly.", default=False, action='store_true')

args = ap.parse_args()

alignFile = open(sys.argv[1])

prevSeqName = ""

def MakeForward(a,b,l, s):
    if (s == "0"):
        return (a,b)
    else:
        return(l-b,l-a)
        
printed = False    
for line in alignFile:
    vals = line.split()
    # parse m1 format   
    #chrX:45799535-45841726|ctg7180000000001/0_68715 chr22 0 1 -3779 81.7782 30466124 30467159 51304566 54761 55791 68715 245
    seqName = vals[0].split("/")[0]
    if (seqName != prevSeqName):
        if (prevSeqName != "" and printed == False):
            sys.stderr.write("Did not find rc for " + prevSeqName + "\n")
        optChrom = vals[1]
        optLen   = int(vals[8])
        optStrand = vals[3]
        (optStart,optEnd) = MakeForward(int(vals[6]), int(vals[7]), optLen, optStrand)
        optLine = line
        if (args.asmgen):
            print optChrom + ":" + str(optStart) + "-" + str(optEnd)
        printed = False
    else:
        curChrom = vals[1]
        curLen   = int(vals[8])
        curStrand = vals[3]
        (curStart,curEnd) = MakeForward(int(vals[6]), int(vals[7]), curLen, curStrand)
        if (curChrom != optChrom):
            continue
        # must be opposite strands
        if (curStrand == optStrand):
            continue

        if (printed == False and (curStart < optStart or curEnd > optEnd)):
            continue

        if (printed == False):
            if (args.length):
                print optChrom + "\t" + str(optStart) + "\t" + str(optEnd) + "\t" + str(curEnd - curStart)
            if (args.gen):
                print curChrom + ":" + str(curStart) + "-" + str(curEnd)
            if (args.asm):
                print seqName + ":" + vals[9] + "-" + vals[10]
            if (args.first):
                print seqName + ":" + vals[9] + "-" + str(int((int(vals[10]) + int(vals[9]))/2))
            if (args.second):
                print seqName + ":" + str(int((int(vals[10]) + int(vals[9]))/2)) + "-" + vals[10]
            printed = True
            
    prevSeqName = seqName


        

    
