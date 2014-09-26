#!/usr/bin/env python

import Tools
import sys
import argparse
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq

ap = argparse.ArgumentParser(description="Print alignments (much) shorter than the expected read length.")
ap.add_argument("sam", help="Input SAM file (s)", nargs="+")
ap.add_argument("--minTrunc", help="Minimal truncated length.", type=int, default=500)
ap.add_argument("--nextHit", help="Write the sequecne of the next hit to this file.", default=None)
ap.add_argument("--minqv", help="Minimum mapping qualiyt", type=int, default=20)
ap.add_argument("--out", help="Write output here.", default=None)

args = ap.parse_args()


starti = 5
endi   = 6
flagi  = 1
tposi  = 3
qnamei = 0
tnamei = 2
seqi   = 8
prevSamEntry = None
prevWasShort = False
references = []
if (args.nextHit is not None):
    nextHitOut = open(args.nextHit, 'w')

if (args.out is not None):
    out = open(args.out, 'w')
else:
    out = sys.stdout


for samFileName in args.sam:
    samFile = open(samFileName, 'r')
    sys.stderr.write(samFileName + "\n")
    for line in samFile:
        if (len(line) > 0 and line[0] == '@'):
            if (line[0:3] == "@SQ"):
                vals = line.split()
                if (vals[1][1:3] == "SN:"):
                    references.append(vals[1][3:])
            continue
        samEntry = Tools.ParseSamLine(line)
        (movie,index,coords) = Tools.ParseReadTitle(samEntry[0])
        mapLen   = samEntry[endi] - samEntry[starti]
        readLen  = coords[1] - coords[0]
        
        if (prevWasShort == True and prevSamEntry[0] == samEntry[0] ):
            curSpan = samEntry[endi] - samEntry[starti]
            differenceBetweenStarts    = prevSamEntry[tposi] - samEntry[tposi]
            differenceBetweenSegments = Tools.GapBetweenIntervals((prevSamEntry[starti], prevSamEntry[endi]),
                                                                   (samEntry[starti], samEntry[endi]))
            if ((prevSamEntry[tnamei] != samEntry[tnamei] or
                 abs(differenceBetweenStarts > curSpan)) and
                 not (samEntry[starti] >= prevSamEntry[starti] and
                      samEntry[starti] <= prevSamEntry[endi] and
                      samEntry[endi] >= prevSamEntry[starti] and
                      samEntry[endi] <= prevSamEntry[endi])):

                overlap = Tools.Overlap((prevSamEntry[starti], prevSamEntry[endi]),
                                        (samEntry[starti], samEntry[endi]))
                curSpan = samEntry[endi] - samEntry[starti]
                prevSpan = prevSamEntry[endi] - prevSamEntry[starti]
                overlapping = False
                if (overlap != 0 and
                    (curSpan / overlap > 0.5 or
                     prevSpan / overlap > 0.5)):
                    overlapping = True
                if (overlapping == False):
                    if (args.nextHit is not None):
                        if (abs(differenceBetweenSegments) < 200):
                            seqname = samEntry[tnamei]+"/" + str(samEntry[tposi]) + "/" + str(len(samEntry[seqi])) + "/" + samEntry[qnamei]
                            seqRecord=Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(samEntry[seqi]), id=seqname,name="",description="")
                            Bio.SeqIO.write(seqRecord, nextHitOut, "fasta")
                    out.write(prevLine)
                    out.write(line)
                
        
        if (readLen - mapLen > args.minTrunc):
            prevWasShort = True
        else:
            prevWasShort = False
        prevSamEntry = samEntry
        prevLine = line

if (args.nextHit is not None):
    nextHitOut.close()

if (args.out is not None):
    out.close()
