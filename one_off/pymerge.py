#!/usr/bin/env python
import sys
prevChr = ""
prevStart = 0
prevEnd   = 0
prevText = ""
endOnPrint = False
inFile = open(sys.argv[1])
for line in inFile:
    vals = line.split()
    if (prevChr == ""):
        (prevChr, prevStart, prevEnd, prevLength, prevText) = (vals[0], int(vals[1]), int(vals[2]), int(vals[3]), "\t".join(vals[4:]))
        continue
    else:
        (chrom, start, end, length, text) = (vals[0], int(vals[1]), int(vals[2]), int(vals[3]), "\t".join(vals[4:]))
        #
        # Case 1, the two lines should not be merged.
        #
        if (chrom != prevChr or start >= prevEnd + 40):
            print "{}\t{}\t{}\t{}\t{}".format(prevChr, prevStart, prevEnd, prevLength, prevText)
            # reset
            endOnPrint = True
            (prevChr, prevStart, prevEnd, prevLength, prevText) = (vals[0], int(vals[1]), int(vals[2]), int(vals[3]), "\t".join(vals[4:]))
        else:
            #
            # case 2, the lines might need to be merged.
            if (chrom == prevChr):
                if (start <= prevEnd + 40):
                    prevEnd = max(prevEnd, end)
                    prevLength = max(length, prevLength)
            endOnPrint = False


if (endOnPrint == True):
    print "{}\t{}\t{}\t{}\t{}".format(prevChr, prevStart, prevEnd, prevLength, prevText)
