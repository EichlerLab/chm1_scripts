#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="Print inversions from bed file created from SamToBed")

ap.add_argument("bed", help="Input bed file")
ap.add_argument("out", help="Output file")
ap.add_argument("--delta", help="Maximum delta of inverted sequence.", type=int, default=0)
ap.add_argument("--minDelta", help="Minimum delta  (to screen for inv+dup)", type=int, default=None)

args=ap.parse_args()
samBed = open(args.bed)
out = open(args.out, 'w')
# parse lines in the format
#(1)chrom (2)start (3)end (4)rname (5)strand (6)readstart (7)readend (8)readlen
#chr15   60754519        60754519        m131125_043202_42215_c100603842550000001823113706171496_s1_p0/108988/11251_19264        1       14571   18769   8013

class Aln:
    def __init__(self, line):
        self.line = line
        vals = line.split()
        if (len(vals) == 8):
            self.chrom = vals[0]
            self.tStart = int(vals[1])
            self.tEnd   = int(vals[2])
            self.rName  = vals[3]
            self.rStrand = int(vals[4])
            self.rStart = int(vals[5])
            self.rEnd   = int(vals[6])
            self.rLen   = int(vals[7])
        else:
            return None



def ComputeDelta(alnA, alnB):
    # First look to see if b is inside of A:
    if (alnB.tStart > alnA.tStart and alnB.tEnd < alnA.tEnd):
        d = alnA.tStart - alnB.tStart + alnB.tEnd - alnA.tEnd
        return d
    
    # Next look to see if there is overlap:
    delta = max(max(alnA.tStart - alnB.tEnd,0), max(alnB.tStart - alnA.tStart, 0))
    return delta


prevAln = None

nProcessed = 0
for line in samBed:
    curAln = Aln(line)
    if (curAln == None):
        continue
    if (prevAln is not None and prevAln.rName == curAln.rName):
        # check for overlap
        if (prevAln.chrom == curAln.chrom):
            delta = ComputeDelta(prevAln, curAln)
            if (delta <= args.delta and ( args.minDelta == None or delta >= args.minDelta) and prevAln.rStrand != curAln.rStrand):
#                out.write(prevAln.line)
                out.write(curAln.line)
    prevAln = curAln
    nProcessed += 1
