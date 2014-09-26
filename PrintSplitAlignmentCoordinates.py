#!/usr/bin/env python

import sys
import argparse
import Tools

ap = argparse.ArgumentParser(description="Print the coordinates of alignments that are split in the genome.")
ap.add_argument("out", help="Output file base  out.ll out.lr, out.rl. out.rr will be created")
ap.add_argument("--sam", help="Input sam file", nargs="+")
ap.add_argument("--allowInter", help="Allow inter chromosomal alignments.", action='store_true', default=False)
ap.add_argument("--allowchrUn", help="Allow inter chromosomal alignments to chrUn.", action='store_true', default=False)
ap.add_argument("--start", help="First sam file to process in the list.", default=0, type=int)
ap.add_argument("--stride", help="Skip this many files.", default=1, type=int)
args = ap.parse_args()

llout = open(outFile + ".ll", "w")
lrout = open(outFile + ".lr", "w")
rlout = open(outFile + ".rl", "w")
rrout = open(outFile + ".rr", "w")


def GetKV(key, val):
    lk = len(key)
    for v in vals:
        if (len(v) >= len(key) and v[0:lk] == key):
            return v[lk:]
    else:
        return None

def GetStrand(value):
    if (value & 16 != 0):
        return 1
    else:
        return 0

#for samFileName in args.sam:

left = 1
right = 2

def SetTrunc(s, e, l):
    if (s > l - e):
        return (left, s)
    else:
        return (right, l-e)

for i in range(args.start, len(args.sam), args.stride):
    samFileName = args.sam[i]
    samFile = open(samFileName)
    prevRead = None
	
    prevLine = ""
    for line in samFile:
        if (len(line) > 0 and line[0] == '@'):
            continue
        vals = line.split()
        curRead = vals[0]
        curChrom = vals[2]
        curPos = int(vals[3])
        curStart = int(GetKV("XS:i:", vals))
        curEnd = int(GetKV("XE:i:", vals))
        curLen = int(GetKV("XL:i:", vals))
        curStrand = GetStrand(int(vals[1]))
        curScore = int(GetKV("AS:i:", vals))
        if (prevRead == curRead):
            prevSam = Tools.SAMEntry(prevLine)
            curSam  = Tools.SAMEntry(line)
            targetOverlap = Tools.Overlap((prevSam.tStart, prevSam.tEnd), (curSam.tStart, curSam.tEnd))
            
            if (prevSam.mapqv > 10 and curSam.mapqv > 10 and
                targetOverlap / (prevSam.tEnd - prevSam.tStart) < 0.1 and
	            targetOverlap / (curSam.tEnd - curSam.tStart) < 0.1):
                if (args.allowInter or curChrom == prevChrom or (args.allowchrUn and (curChrom[0:5] == 'chrUn' or prevChrom[0:5] == 'chrUn'))):
                    c = (curChrom, curStrand, curSam.tStart, curSam.tEnd, curStart, curEnd, curLen)
                    p = (prevChrom, prevStrand, prevSam.tStart, prevSam.tEnd, prevStart, prevEnd, prevLen)
                    if (curSam.tStart > prevSam.tStart):
                        t = c
                        c = p
                        p = t

                    
                    # First determine which end is truncated on the first fragment (p)
                    (prevTrunc, prevTruncLen) = SetTrunc(p[4], p[5], p[6])
                    (curTrunc, curTruncLen)   = SetTrunc(c[4], c[5], c[6])

                    
                    
                    outFile.write( curRead + "\t" + str(curScore) + "\t" + str(c[1]) + "\t" + c[0] + "\t" +  str(c[2]) + "\t" + str(c[3]) + "\t" + str(prevScore)+ "\t" + str(p[1]) + "\t" + str(p[0]) + "\t" + str(p[2]) + "\t" + str(p[3]) + "\n")
	#            print str((curStart, prevStart, curEnd, prevEnd))
	#            print "disjoint alignments: "
	#            print line
	#            print prevLine
        prevRead = curRead
        prevChrom = curChrom
        prevStart = curStart
        prevEnd   = curEnd
        prevLine = line
        prevStrand = curStrand
        prevLen = curLen
        prevPos = curPos
        prevScore = curScore
