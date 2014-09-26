#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import DNAAlphabet

import argparse
ap = argparse.ArgumentParser(description="Add an additional column to a counts file describing if a sequence is in a flanking or region of interest.")
ap.add_argument("counts", help="Counts file.")
ap.add_argument("seq", help="Sequence file. This should be referenced by the first column of the counts file.")
ap.add_argument("flank", help="How much from the side of each sequence is flanking", type=int)
ap.add_argument("output", help="Output file.")
args = ap.parse_args()

countsFile = open(args.counts, 'r')
outFile = open(args.output, 'w')

seqFile = open(args.seq, 'r')
seqDict = SeqIO.to_dict(SeqIO.parse(seqFile, "fasta"))

for key,seq in seqDict.iteritems():
    seqDict[key] = seq.upper()

for line in countsFile:
    vals = line.split()
    seqName = vals[0]
    state = 0
    if (seqName not in seqDict):
        state = 2
    else:
        kmer = vals[1]
        kmerrc = Seq.Seq(kmer, DNAAlphabet()).reverse_complement()
        kmerPos = seqDict[seqName].seq.find(kmer)
        kmerrcPos = seqDict[seqName].seq.find(kmerrc)
        pos = kmerPos
        if (kmerrcPos >= 0):
            pos = kmerrcPos
        if (pos >= 0):
            distToFlank = min(pos, len(seqDict[seqName].seq) - pos)
#            print str(pos) + " " + str(distToFlank) + " " + str(len(seqDict[seqName].seq)) + " " + seqName
            if (distToFlank < args.flank):
                state = 1
            else:
                state = 0
        else:
            state = 2
    outFile.write(line.strip() + "\t" + str(state) + "\n")

outFile.close()

