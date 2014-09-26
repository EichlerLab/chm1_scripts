#!/usr/bin/env python
import sys
prevChr = ""
prevStart = 0
prevEnd   = 0
prevText = ""
endOnPrint = False

for line in sys.stdin.readlines():
    vals = line.split()
    if (prevChr == ""):
        (prevChr, prevStart, prevEnd, prevRepeat, prevExpand) = (vals[0], int(vals[1]), int(vals[2]), vals[3], int(vals[8]))
#        (prevChr, prevStart, prevEnd, prevLine, prevExpand) = (vals[0], int(vals[1]), int(vals[2]), "\t".join(vals[0:7]), int(vals[10]))
        continue
    else:
        (chrom, start, end, repeat, expand) = (vals[0], int(vals[1]), int(vals[2]), vals[3], int(vals[8]))
        if (chrom != prevChr or start != prevStart):
            print "{}\t{}\t{}\t{}\t{}\t{}".format(prevChr, prevStart, prevEnd, prevRepeat, int(prevEnd) - int(prevStart), prevExpand)
            # reset
            endOnPrint = True
            (prevChr, prevStart, prevEnd, prevRepeat, prevExpand) = (vals[0], int(vals[1]), int(vals[2]), vals[3], int(vals[8]))
        else:
            if (chrom == prevChr):
                prevExpand += int(vals[8])
            endOnPrint = False


if (endOnPrint == True):
    print "{}\t{}\t{}\t{}\t{}\t{}".format(prevChr, prevStart, prevEnd, prevRepeat, int(prevEnd) - int(prevStart), prevExpand)
