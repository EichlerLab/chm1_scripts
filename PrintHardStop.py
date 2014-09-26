#!/usr/bin/env python
import sys
import argparse
import Tools
import os

ap = argparse.ArgumentParser(description="Print hard-stop coordinates, with a little wiggle-room, and the read names.")
ap.add_argument("sam", help="Input sam files.", nargs="+")
ap.add_argument("--out", help="Output file, 'stdout' defaults to stdout", default="stdout")
ap.add_argument("--window", help="Print this window around a hard-stop sequence.", default=100)
ap.add_argument("--delta", help="Minimum truncation to trigger a hard-stop", default=500)
ap.add_argument("--minAlign", help="Minimum alignment length to consider.", default=800)
ap.add_argument("--splitread", help="Write splitread alignments to this file.", default=None)
args = ap.parse_args()
titlei = 0
flagi = 1
tnamei = 2
tposi  = 3
starti = 5
endi = 6
mapqvi = 4
tleni = 9

if (args.out == "stdout"):
    outFile = sys.stdout
else:
    outFile = open(args.out, "w")

if (args.splitread is not None):
    if (args.out == "stdout"):
        print "ERROR, cannot have splitread alignment file with stdout."
    splitReadFileName = args.splitread
    splitReadFile = open(splitReadFileName, 'w')
    
index = 0
prevEntry = None
nhs = 0
nlines = 0

left = 0
right = 1

def IsHardStop(line, minClip):
    entry = Tools.ParseSamLine(line)
    cigar =  line.split()[5]
    ops, lengths = Tools.CIGARToArrays(cigar)
    (preSoftClip, sufSoftClip) = Tools.GetSoftClip(ops, lengths)
#    if (Tools.GetStrand(entry[flagi]) == 1):
#        # swap entries.
#        (preSoftClip, sufSoftClip) = (sufSoftClip, preSoftClip)

    alnStart = entry[starti] + preSoftClip
    alnEnd   = entry[endi]   - sufSoftClip
    if (rEnd - rStart < args.minAlign):
        return None

    # If both sides are clipped, that is hard to interpret, and so dkip
    if (preSoftClip > minClip and sufSoftClip > minClip):
        return None

    if (preSoftClip > minClip):
        return [entry[tnamei], entry[tposi], entry[tposi]+ entry[tleni], alnStart, alnEnd, preSoftClip, entry[titlei], "left", Tools.GetStrand(entry[flagi])]
    if (sufSoftClip > minClip):
        return [entry[tnamei], entry[tposi], entry[tposi]+ entry[tleni], alnStart, alnEnd, sufSoftClip, entry[titlei], "right", Tools.GetStrand(entry[flagi])]


nhs += 1
curTitle = ""
prevTitle = ""
import pdb
for samFileName in args.sam:
    if (os.path.exists(samFileName) == False):
        continue
    samFile = open(samFileName, 'r')
    index +=1
    lines = []
    for line in samFile:
        if (line[0] == "@"):
            continue
        else:
            entry = Tools.ParseSamLine(line)
            order = "primary"

            if (entry[mapqvi] < 30):
                continue

            if (entry[flagi] & 256 != 0):
                order = "secondary"
            nlines += 1
            (barcode, zmw, (rStart, rEnd)) = Tools.ParseReadTitle(entry[titlei])

            curTitle = entry[titlei]

            if (curTitle != prevTitle):
                hs = [ IsHardStop(l, args.delta) for l in lines ]
                numHS = len(hs) - hs.count(None)
                if (numHS == 1 or numHS == 2 and args.splitread is None):
#                    pdb.set_trace()
#                    print str(hs)
                    for h in hs:
                        if ( h is not None):
                            outFile.write("\t".join([str(i) for i in h]) + "\n")
                elif (numHS == 2 and args.splitread is not None):
                    splitReadFile.write("\t".join([str(i) for i in hs[0]]) + "\t" + "\t".join([str(i) for i in hs[1]])+"\n")
                lines = []
            prevTitle = curTitle
            lines.append(line)
            prevEntry = entry
            if (nlines % 10000 == 0):
                sys.stderr.write("{} files {} lines {} hs \n".format(index, nlines, nhs))
if (outFile != sys.stdout):
    outFile.close()
                            
                
