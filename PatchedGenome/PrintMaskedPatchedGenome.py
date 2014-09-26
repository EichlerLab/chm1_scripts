#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord


ap = argparse.ArgumentParser()
ap.add_argument("genome", help="Reference genome.")
ap.add_argument("patches", help="Inserted sequences.")
ap.add_argument("outputGenome", help="Output genome")
ap.add_argument("--flank", help="Add this extra flanking sequence.", default=100, type=int)
ap.add_argument("--chrom", help="Patch this chromosome only.", default=None)

args = ap.parse_args()

# assume the patches file is sorted by chromosome
patchesByChrom = {}

patchesFile = open(args.patches)

for line in patchesFile:
    vals = line.split()
    if (vals[0] not in patchesByChrom):
        patchesByChrom[vals[0]] = []
    patchesByChrom[vals[0]].append([max(0, int(vals[1])-args.flank), int(vals[2])+args.flank, vals[3], int(vals[1]), int(vals[2])])



genome = open(args.genome)
modGenome = open(args.outputGenome, 'w')
patches = []
patchNames = []

for chrom in SeqIO.parse(genome, "fasta"):
    #
    # Due to python's string immutability, in order to quickly mask,
    # create a list of substrings that will be joined with
    # interleaving spans of N's.
    # This requires non-overlapping windows.  Merge them first.
    #
    if (args.chrom is not None and args.chrom != chrom.id):
        continue

    maskOut = [(0,0)]
    name = chrom.id
    if (name in patchesByChrom):
        curStart = patchesByChrom[name][0][0]
        curEnd   = patchesByChrom[name][0][1]

        for i in range(0,len(patchesByChrom[name])-1):
            if (patchesByChrom[name][i][1] >=
                patchesByChrom[name][i+1][1]):
                curEnd = patchesByChrom[name][i+1][1]
            else:
                maskOut.append((curStart, curEnd))
                curStart = patchesByChrom[name][i+1][0]
                curEnd   = patchesByChrom[name][i+1][1]
        lasti = len(patchesByChrom[name])-1

        if (curStart == patchesByChrom[name][lasti][0]):
            maskOut.append((curStart, curEnd))

        maskOut.append((len(chrom.seq), len(chrom.seq)))

        # Prior to masking the genome, create the list of insertions.
        chromSeq = chrom.seq.tostring()
        patches = patches + [ chromSeq[patchesByChrom[name][i][0]:patchesByChrom[name][i][3]] + patchesByChrom[name][i][2] + chromSeq[patchesByChrom[name][i][1]:patchesByChrom[name][i][4]] for i in range(0,len(patchesByChrom[name]))]
        patchNames = patchNames + ["{}:{}-{}/patch".format(name, patchesByChrom[name][i][3], patchesByChrom[name][i][4]) for i in range(0,len(patchesByChrom[name]))] 

    
        # now maskOut should be non-overlapping coordinates.
        maskedContig = ''.join([ chromSeq[maskOut[i-1][1]:maskOut[i][0]] + 'N'*(maskOut[i][1] - maskOut[i][0]) for i in range(1,len(maskOut))])
    
        SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(maskedContig), id=chrom.id, name="", description=""), modGenome, "fasta")
    else:
        SeqIO.write(chrom, modGenome, "fasta")

for i in range(0,len(patches)):
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(patches[i]), id=patchNames[i], name="", description=""), modGenome, "fasta")

modGenome.close()
genome.close()
