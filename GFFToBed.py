#!/usr/bin/env python
import sys

if (len(sys.argv) < 3):
    print "usage: GFFToBed.py in.gff out.bed.  Prints insertions or deletions."
    sys.exit(0)
    
inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
for line in inFile:
    vals = line.split()
        
    if (len(vals) == 9 and (vals[2] == "insertion" or vals[2] == "deletion")):
        chrom = vals[0]
        start = vals[3]
        end   = vals[4]
        varSeq = vals[8].split(';')[1].split('=')[1]
        varLen = max(len(varSeq), int(end) - int(start))
        outFile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, vals[2], varLen, varSeq))

outFile.close()
