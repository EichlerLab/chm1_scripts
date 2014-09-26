#!/usr/bin/env python

import Tools

import argparse

ap = argparse.ArgumentParser(description="Print aligned blocks in a sam file.")
ap.add_argument("sam", help="sam file.")
ap.add_argument("block", help="Block file, in bed format.")
ap.add_argument("--maxgap", help="Maximum gap to merge in a block.", default=20)

args = ap.parse_args()

outFile = open(args.block, 'w')

sf = open(args.sam)
index = 0
for line in sf:
    if (line[0] == '@'):
        continue
    else:
        sam = Tools.SAMEntry(line)
        qPos = sam.qStart
        tPos = sam.tStart
        i = 0
        nCigar = len(sam.lengths)
        while (i < nCigar and (sam.ops == 'S' or sam.ops == 'H')):
            i += 1
        
        while (i < nCigar):
            # Make sure every block starts on a match
            while (i < nCigar and (sam.ops[i] == 'I' or sam.ops[i] == 'D')):
                if (sam.ops[i] == 'I'):
                    qPos += sam.lengths[i]
                elif (sam.ops[i] == 'D'):
                    tPos += sam.lengths[i]
                i+=1

            qStart = qPos
            qEnd = qPos
            tEnd = tPos
            nIns = 0
            nDel = 0
            while (i < nCigar and (sam.ops[i] == 'M' or ( (sam.ops[i] == 'I' or sam.ops[i] == 'D') and sam.lengths[i] < args.maxgap))):
                if (sam.ops[i] == 'M'):
                    qEnd += sam.lengths[i]
                    tEnd += sam.lengths[i]
                if (sam.ops[i] == 'I'):
                    qEnd += sam.lengths[i]
                    nIns += sam.lengths[i]
                if (sam.ops[i] == 'D'):
                    tEnd += sam.lengths[i]
                    nDel += sam.lengths[i]
                i += 1
            outFile.write("{}\t{}\t{}\t{:2.2f}\t{}\n".format(sam.title, index, qStart, qEnd, float(abs((qEnd - qStart) - (nIns+nDel)))/(qEnd-qStart)))
            tPos = tEnd
            qPos = qEnd
        index += 1
