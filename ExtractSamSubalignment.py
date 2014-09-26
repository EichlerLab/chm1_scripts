#!/usr/bin/env python


import argparse
import Tools
import sys
import pdb
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq 

ap = argparse.ArgumentParser(description="Given an interval on the reference, extract the portion of the cigar covering that interval.")
ap.add_argument("sam", help="Input SAM file.")
ap.add_argument("interval", help="Interval on reference")
ap.add_argument("--wiggle", help="Extend interval by this amount", type=int, default=0)
ap.add_argument("--out", help="Output fasta file. For now this is just printing fasta.", default="/dev/stdout")
ap.add_argument("--extName", help="Add this tag to the sequence name, helpful for delimiting sequence source in downstream analysis.", default="")
ap.add_argument("--append", help="Append output instead out create new", action='store_true', default=False)

args = ap.parse_args()

if (args.extName != ""):
    args.extName = "/" + args.extName
    
region = Tools.ParseRegionStr(args.interval)
expanded = (region[0], region[1] - args.wiggle, region[2] + args.wiggle)
region = expanded

samFile = open(args.sam)

openMode = 'w'
if (args.append):
    openMode = 'a'

outFile = open(args.out, openMode)

# assume any line in the sam file can be a b

samHeader = Tools.SAMHeader()

line = samHeader.ReadHeader(samFile)

while (line != ""):
    samEntry = Tools.SAMEntry(line)
    if (samEntry.tName == region[0]):
        if (not ((samEntry.tStart <= region[1] and samEntry.tEnd >= region[1]) or
                 (samEntry.tStart <= region[2] and samEntry.tEnd >= region[2]) or
                 (samEntry.tStart <= region[1] and samEntry.tEnd >= region[2]))):
            line = samFile.readline()
            continue
        else:
            pass
    else:
        line = samFile.readline()
        continue
    
    #
    # The alignment overlaps the interval.
    #

    # advance until overlapping with qpos
    #
    i = 0

    while (i < len(samEntry.ops) and (samEntry.ops[i] == 'S' or samEntry.ops[i] == 'H')):
        i += 1
    e = len(samEntry.ops)
    while (e > 0 and (samEntry.ops[e-1] == 'S' or samEntry.ops[e-1] == 'H')):
        e -= 1
    ops = samEntry.ops[i:e]
    lengths = samEntry.lengths[i:e]
    
    tPos = [samEntry.tStart]*(len(ops)+1)
    qPos = [samEntry.qStart]*(len(ops)+1)
    tOp  = ['M']*(len(ops)+1)
    totalM = 0
    for i in range(0,len(ops)):
        if (ops[i] == 'M'):
            tPos[i+1] = tPos[i] + lengths[i]
            qPos[i+1] = qPos[i] + lengths[i]
            totalM += lengths[i]
        elif (ops[i] == 'D'):
            tPos[i+1] = tPos[i] + lengths[i]
            qPos[i+1] = qPos[i]
            tOp[i] = 'D'
        elif (ops[i] == 'I'):
            qPos[i+1] = qPos[i] + lengths[i]
            tPos[i+1] = tPos[i]
    # Now find the interval in the t pos.
    s = 0
    while (s < len(tPos)-1 and (not (region[1] >= tPos[s] and region[1] < tPos[s+1]))):
        s += 1
    e = s
    while (e < len(tPos)-1 and (not (region[2] >= tPos[e] and region[2] < tPos[e+1]))):
        e += 1
        

    # The interval should be from s to e
    if (s == e and e == len(tPos)):
        print "ERROR. Did not find an interval that should have been in tStart...tEnd."
        sys.exit(0)

        startRegionOffset = max(0,region[1] - tPos[s])
#    pdb.set_trace()    
    prefixDeletion = 0
    if (tOp[s] == 'D'):
        prefixDeletion = lengths[s]
    
    endRegionOffset   = region[2] - tPos[e]
    if (tOp[e] == 'D'):
        endRegionOffset -= lengths[s]

    sys.stderr.write("tpos starting around " + str(tPos[s]) + "\n")
    sys.stderr.write("qpos starting around " + str(qPos[s]) + "\n")
    sys.stderr.write("offset " + str(startRegionOffset) + " total q: " + str(qPos[s] + startRegionOffset) + "\n")
    seqPos = qPos[s] + startRegionOffset
    queryGapSum = sum(lengths[s:e])
    regionLength = endRegionOffset + queryGapSum - startRegionOffset - prefixDeletion
    
#    pdb.set_trace()
#    print "got interval " + str(s) + "\t" + str(e) + "\t" + str(regionLength) +"\t" +str(region[2] - region[1]) + "\t" + str(region[2] - region[1] - regionLength)
    seq = samEntry.seq[seqPos:seqPos+regionLength]
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq), id="{}:{}-{}{}".format(region[0],region[1],region[2], args.extName), name="", description=""), outFile, "fasta")
    line = samFile.readline()


outFile.close()
