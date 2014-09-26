#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="""Take an inversion duplication table in the format
insertion chr1:6043687-6086793/22279/22399/0_120 ctg7180000000003_quiver 0 0 -600 100 11880 12000 31983 0 120 120 119
insertion chr1:6043687-6086793/22279/22399/0_120 ctg7180000000003_quiver 0 0 -589 99.1667 12748 12868 31983 0 120 120 95
(CACCC)n,108,17.4
and turn it into calls dup_len dup_pos""")

ap.add_argument("table", help="Candidate invdup table.")
ap.add_argument("--plotCommands", help="Write plotcommands here.")
ap.add_argument("--minIdentity", help="Minimum percent identity (default=98.0)", default=98.0,type=float)
args=ap.parse_args()

tabFile = open(args.table)

def ValsToForward(v):
    return ToForward(v[4],v[7],v[8],v[9])

def ToForward(strand,s,e,l):
    if (strand == "0"):
        return (int(s),int(e), int(l))
    else:
        return (int(l)-int(e),int(l)-int(s), int(l))

line = tabFile.readline()
if (line != '' and line.strip() != "START"):
    print "malformatted invdup table file"
    sys.exit(1)

plotCommandFile = None    
if (args.plotCommands is not None):
    plotCommandFile = open(args.plotCommands, 'w')

while (line != '' and tabFile):
    assemblyFile = tabFile.readline()
    regionLine   = tabFile.readline()
    repeats = []
    while (True):
        line = tabFile.readline()
        line = line.strip()
        if (line == '' or line == "START"):
            break
        alnLine2 = tabFile.readline().strip()
        repeatLine = tabFile.readline().strip()
        repeats.append( (line, alnLine2, line, repeatLine) )
        # ends when line == START
        
    if (len(repeats) == 1):
        v1 = repeats[0][0].split()
        v2 = repeats[0][1].split()
        ident1 = float(v1[6])
        ident2 = float(v2[6])
        if (ident1 < args.minIdentity or ident2 < args.minIdentity):
            continue
        region = regionLine.strip()
        print "{} {} {:2.2f} {:2.2f}".format(assemblyFile.strip(), region, ident1, ident2)
        r1c = ValsToForward(v1)
        r2c = ValsToForward(v2)
        start = min(r1c[0], r2c[0])
        end   = max(r1c[1], r2c[1])

        rgnChr = region.split(":")[0]
        [rgnStart, rgnEnd] = [ int(i) for i in region.split(":")[1].split("-")]
        if (plotCommandFile is not None):
           plotCommandFile.write("$PBS/DotPlot.py --query {} --target /var/tmp/mchaisso/ucsc.hg19.fasta --region {}:{}-{} --queryStart {} --queryEnd {} --savefig {}_{}-{}.png --nolegend --matches dot:13 &\n".format(assemblyFile.strip(), rgnChr, rgnStart + start - 5000, rgnEnd - 5000, start - 10000, end + 10000,  rgnChr, rgnStart, rgnEnd))


if (plotCommandFile is not None):
    plotCommandFile.close()
        
