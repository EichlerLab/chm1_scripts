#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="print overlapping alignments.")
ap.add_argument("input", help="input paired read bed file. adjacent lines are from the same read.")
ap.add_argument("--window", help="Maximum allowed window", default=20000, type=int)
ap.add_argument("--unique_subreads", help="Only print support from one subread.", action='store_true', default=False)
ap.add_argument("--mingap", help="This many bases on the read must not be aligned.", default=100, type=int)

args = ap.parse_args()

inFile = open(args.input)

reads={}
import pdb
nTGap = 0
while True:
    l1 = inFile.readline()
    l2 = inFile.readline()
    if (l1 == "" or l2 == ""):
        break
    v1 = l1.split()
    v2 = l2.split()
    intv1 = [v1[0], int(v1[1]), int(v1[2])]
    intv2 = [v2[0], int(v2[1]), int(v2[2])]

    
    overlapStart = max(intv1[1], intv2[1])
    overlapEnd   = min(intv1[2], intv2[2])
    overlap=max(0, overlapEnd - overlapStart)
    overlapPercent = float(overlap) / (intv2[2] - intv2[1])

    subreadCoordinates = v1[3].split("/")[-1].split("_")
    subreadLength = int(subreadCoordinates[1]) - int(subreadCoordinates[0])

    v1ts = int(v1[1])
    v1te = int(v1[2])
    v2ts = int(v2[1])
    v2te = int(v2[2])
    
    v1qs = int(v1[6])
    v1qe = subreadLength - int(v1[7])
    v2qs = int(v2[6])
    v2qe = subreadLength - int(v2[7])

    doPrint = False
    #
    # Make sure the target alignments aren't too far away.
    #
    tGap = max(v2ts - v1te, v1ts - v2te)
    qGap = max(v2qs - v1qe, v1qs - v2qe)

    #
    # Determine the breakpoints
    #
    if (tGap < 0):
        continue
    
    if (tGap == v2ts - v1te):
        tBP = (v1te, v2ts)
    else:
        tBP = (v2te, v1ts)



    if (tBP[0] > tBP[1]):
        pdb.set_trace()

    if (tGap < args.window):
        doPrint = True

        
    #
    # Make sure the query alignments aren't too close.
    #

    if (qGap < args.mingap):
        doPrint = False

    if (qGap == v2qs - v1qe):
        qBP = (v1qe, v2qs)
    else:
        qBP = (v2qe, v1qs)

    if (args.unique_subreads):
        read = "/".join(v1[3].split('/')[0:2])
        if read in reads:
            doPrint = False
        else:
            reads[read] = True
    #
    # The two alignments must be in the same orientation (so that the gap is inverted)
    #
    if (v1[8] != v2[8]):
        doPrint = False
        
    if (doPrint):
        #
        # Print both alignments for now
        #
        l1 = l1.strip() + "\t{}\t{}\n".format(qGap, tGap)
        l2 = l2.strip() + "\t{}\t{}\n".format(qGap, tGap)
        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(v1[0], tBP[0], tBP[1], v1[3], qBP[0], qBP[1], tBP[1]-tBP[0], qBP[1] - qBP[0]))
