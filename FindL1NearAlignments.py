#!/usr/bin/env python

import sys
import argparse
import Tools
import os

ap = argparse.ArgumentParser(description="Find L1s near alignments.")
ap.add_argument("l1", help="L1 alignments in m1 format.")
ap.add_argument("sam", help="Input sam files.", nargs="+")
ap.add_argument("--window", help="Print this window around sides with line.",type=int, default=100)
ap.add_argument("--delta", help="Restrict line matches to be this close to alignment ends.", type=int, default=500)
ap.add_argument("--covered", help="If this l1 is fully contained in an alignment, print to a file.", default=None)
ap.add_argument("--notcovered", help="If this l1 is not covered by an alignment, print to this file.", default=None)

args = ap.parse_args()


l1File = open(args.l1)
l1Alignments = {}
if (args.covered is not None):
    coveredFile = open(args.covered, 'w')
if (args.notcovered is not None):
    notCoveredFile = open(args.notcovered, 'w')
    
for line in l1File:
    vals = line.split()
    if (int(vals[3]) == 1):
        qlen = int(vals[11])
        qs   = qlen - int(vals[10])
        qe   = qlen - int(vals[9])
        vals[9] = qs
        vals[10] = qe
    l1Alignments[vals[0]] = vals[1:]


for samFileName in args.sam:
    if (os.path.exists(samFileName) == False):
        continue
    
    samFile = open(samFileName)
    sys.stderr.write(samFileName + "\n")
    for line in samFile:
        if (line[0] == '@'):
            continue
        samVals = Tools.ParseSamLine(line)
        if (samVals[Tools.mapqvi] < 30):
            continue
        if (samVals[Tools.titlei] in l1Alignments):
            l1Align = l1Alignments[samVals[Tools.titlei]]
            # m1 format:
            #0                                                                       1    2 3 4     5       6ts7te 8tl  9qs  10qe 11qlen  12clust-score
            #m130928_232712_42213_c100518541910000001823079209281310_s1_p0/59/0_7808 line 0 0 -8642 82.7235 5 2555 5403 5326 7754 7916 1043
            qs = int(l1Align[9])
            qe = int(l1Align[10])
            doPrint = False
            alns = samVals[Tools.qstarti]
            alne = samVals[Tools.qendi]
            if (qs >= alns and qe <= alne):
                if (args.covered is not None):
                    coveredFile.write(line)
            else:
                if (args.notcovered is not None):
                    notCoveredFile.write(line)
                    
#            print str(((qs,qe),(samVals[Tools.qstarti],samVals[Tools.qendi])))
            if (qe > samVals[Tools.qstarti] - args.delta and qe < samVals[Tools.qstarti]):
                windowStart = max(samVals[Tools.tpos] - args.window, 0)
                windowEnd   = samVals[Tools.tpos]
                
                doPrint = True
            if (qs > samVals[Tools.qendi] and qs < samVals[Tools.qendi] + args.delta):
                windowStart = samVals[Tools.tposi] + samVals[Tools.tleni]
                windowEnd   = samVals[Tools.tposi] + samVals[Tools.tleni] + args.window
                doPrint = True
            if (doPrint):
                print "{}\t{}\t{}\t{};{};{};{}".format(samVals[Tools.tnamei], windowStart, windowEnd,
                                                       samVals[Tools.qstarti], samVals[Tools.qendi], qs, qe)
