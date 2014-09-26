#!/usr/bin/env python

import pysam
import sys
import argparse
import Tools
import Align
import pdb
from  Bio import SeqIO


ap = argparse.ArgumentParser(description="Print gaps in a SAM file.")
ap.add_argument("genome", help="Genome file with a .fai")
ap.add_argument("sam", help="Sam file of alignment.")
ap.add_argument("--minLength", help="Minimum gapLength.", default=30, type=int)
ap.add_argument("--outFile", help="Print output here, default= stdout", default=None)
ap.add_argument("--context", help="Print surrounding context", default=0, type=int)
ap.add_argument("--condense", help="Pack indels if the matches separating them is less than this value.", default=0, type=int)
ap.add_argument("--tsd", help="Attempt to find Target Site Duplications at most this length", default=0, type=int)

args = ap.parse_args()

genome = file(args.genome, 'r')
handle = open(args.genome, "r")


#samFile = pysam.Samfile(args.sam)
samFile = open(args.sam)

if (args.outFile is None):
    outFile = sys.stdout
else:
    outFile = open(args.outFile, 'w')


    
genomeDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


fai = Tools.ReadFAIFile(args.genome + ".fai")

M = 'M'
I = 'I'
D = 'D'
N = 'N'
S = 'S'
H = 'H'
P = 'P'

    
for line in samFile.readlines():
    if (line[0] == "@"):
        continue

    aln = Tools.SAMEntry(line)



    if (aln.flag & 256 != 0):
        continue
    if (aln.mapqv < 10):
        continue

    tPos = aln.tStart
    qPos = 0
    #
    # condense matches.
    #
    packedCigar = []
    i = 0
    i1 = 1
    niter = 0
    maxGap = 0
    maxGapType = 0
    packedOps = []
    packedLengths = []
#    print str(aln.ops)
    while (i < len(aln.lengths)):
        l = aln.lengths[i]
        op = aln.ops[i]
#        print str((l,op))
#        c = aln.cigar[i]
        j = i
#        sys.stderr.write(str(i) + "\n")
        if (op == I or op == D):
            if (l > maxGap):
                maxGap = l
                maxGapType = op
                
        if (op == I or op == D and i < len(aln.ops) - 2 and aln.ops[i+2][0] == op):
            matchLen = 0
            gapLen   = 0
            while (j+2 < len(aln.ops) and aln.ops[j+2][0] == op and aln.ops[j+1][0] == M and aln.lengths[j+1] < args.condense):

                matchLen += aln.lengths[j+1]
                gapLen   += aln.lengths[j+2]
                j+=2
            if (j > i):
                newIndel = (op, l+gapLen)
                newMatch = (M, matchLen)
                packedOps.append(op)
                packedLengths.append(l+gapLen)

                packedOps.append(M)
                packedLengths.append(matchLen)

#                sys.stderr.write("advance : " + str(j) + "\n")
            else:
                packedLengths.append(l)
                packedOps.append(op)

        else:
            packedLengths.append(l)
            packedOps.append(op)

        i = j + 1
        niter +=1
        if (niter > len(aln.ops)):
            print "ERROR! too many interations."
    sys.stdout.write(aln.title + "\t" + str(maxGap) + "\t" + str(maxGapType) + "\n")
    for i in range(len(packedOps)):
        op = packedOps[i]
        l  = packedLengths[i]

        if (op == M or op == N or op == S):
            tPos += l
            qPos += l
        if (op == I):
            if (l > args.minLength):
                chrName = aln.tName
#                chrName = samFile.references[aln.tid]

                gapSeq = aln.seq[max(0,qPos-args.context):min(qPos+l+args.context, len(aln.seq))]
                tsd = "notsd"
                if (len(gapSeq) == 0):

                    print "ERROR, gap seq is of zero length"
                if (args.tsd):
                    # try and find the target site duplications, this may be on either side of the alignemnt
                    tsdSuffix = gapSeq[-args.tsd:]
                    tsdSuffix = tsdSuffix.upper()
                    tsdPrefix = gapSeq[0:args.tsd]
                    tsdPrefix = tsdPrefix.upper()
                    targetPrefix = genomeDict[chrName].seq[tPos-args.tsd:tPos]
                    targetPrefix = targetPrefix.upper()
                    targetSuffix = genomeDict[chrName].seq[tPos:tPos+args.tsd]
                    targetSuffix = targetSuffix.upper()
                    (sp, ss, sScore) = Align.SWAlign(tsdSuffix, targetPrefix)
                    (pp, ps, pScore) = Align.SWAlign(tsdPrefix, targetSuffix)
#                    print tsdSuffix + " " + tsdPrefix + " " + targetPrefix + " " + targetSuffix
                    if (sScore > pScore ):
                        tsd = ss
                    elif (pScore > sScore ):
                        tsd = ps
                    
                outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrName, tPos, tPos + l, "insertion", l, gapSeq, tsd, aln.title))
            qPos += l
        if (op == D):
            if (l > args.minLength):
#                chrName = samFile.references[aln.tid]
                chrName = aln.tName
                if (tPos > fai[chrName][0]):
                    print "ERROR! tpos is past the genome end." + str(tPos) + " " + str(fai[chrName][0])
                delStart = max(tPos - args.context, 0)
                delEnd   = min(tPos + args.context + l, fai[chrName][0])
                delSeq = genomeDict[chrName].seq[delStart:delEnd].tostring()
#                delSeq = Tools.ExtractSeq([chrName, delStart, delEnd], genome, fai)
                outFile.write("{}\t{}\t{}\t{}\t{}\t{}\tno_tsd\t{}\n".format(chrName, tPos, tPos + l, "deletion", l, delSeq, aln.title))
            tPos += l
        if (op == H):
            pass
#            tPos += l

        
outFile.close()
