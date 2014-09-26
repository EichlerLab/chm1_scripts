#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import pdb

ap = argparse.ArgumentParser(description="Create a new genome from a patch-bed file where the inserted sequences are either spliced in, or inserted (overwriting) existing sequence.")
ap.add_argument("genome", help="Input genome.")
ap.add_argument("patches", help="Patched bed file.  This is in the format:\nchr\tstart\tend\tseq\n   In 'patch' mode, the genome will be deleted from start to end, and 'seq' inserted at 'start'.  In 'insert' mode, the sequence will be pasted over the genome starting at 'start'.")
ap.add_argument("--patchedGenome", help="Output genome with patches. In 'patch' mode, a file 'patchedGenome'.chain will be created that is compatible with liftOver to map annotations from the reference genome to the patched genome.", default=None)
ap.add_argument("--chain", help="Write chain to this file.", default=None)
ap.add_argument("--liftedPatch", help="Write the bed file out with the coordinates in the new file.", default=None)
ap.add_argument("--mode", help="Type of new genome to create: patch|insert.  If a patch file is creatd, ", choices=['patch','insert'], default='patch')


args = ap.parse_args()


patchesFile = open(args.patches)
patchesByChrom = {}

for line in patchesFile:
    vals = line.split()
    if (vals[0] not in patchesByChrom):
        patchesByChrom[vals[0]] = []
    patchesByChrom[vals[0]].append([int(vals[1]), int(vals[2]), vals[3]])
    


genome = open(args.genome)
if (args.patchedGenome is not None):
    modGenome = open(args.patchedGenome, 'w')

if (args.liftedPatch is not None):
    liftBed = open(args.liftedPatch, 'w')
    
chainFile = None
if (args.mode == 'patch'):
    if (args.chain is not None):
        chainFile = open(args.chain, 'w')
    else:
        chainFile = open(args.patchedGenome + ".chain", 'w')
    
patches = []
patchNames = []
chainId = 1

for chrom in SeqIO.parse(genome, "fasta"):

    if (chrom.id not in patchesByChrom):
        #
        # for now skip any chromosomes not modified.  Later when
        # testing is complete they should be output with some sort of
        # no-op for the chain.
        #
        continue

    #
    # Start cutting the genome from pos 0.
    # 
    prevChromEnd = 0
    name = chrom.id

    curStart = patchesByChrom[name][0][0]
    curEnd   = patchesByChrom[name][0][1]
    chromIntervals = []
    chromInsertions = []
    for i in range(0,len(patchesByChrom[name])):
        if (args.mode == 'insert'):
            # insert mode.  Make sure this insertion is not going to overwrite the next one.
            #
            insertionSeq = patchesByChrom[name][i][2]
            if (i < (len(patchesByChrom[name])-1) and patchesByChrom[name][i][0] + len(insertionSeq) > patchesByChrom[name][i+1][0]):
                # the insertion is too long, truncate it.
                overlap =  (patchesByChrom[name][i][0] + len(insertionSeq) ) - patchesByChrom[name][i+1][0]
                insertionSeq = insertionSeq[0:-overlap]
                patchesByChrom[name][i][2] = insertionSeq

            patchesByChrom[name][i][1] = patchesByChrom[name][i][0] + len(insertionSeq)

                
        chromIntervals.append((prevChromEnd, patchesByChrom[name][i][0]))
        chromInsertions.append(patchesByChrom[name][i][2])
        prevChromEnd = patchesByChrom[name][i][1]


    # add a patch from the last insertion to teh end of the genome
    chromIntervals.append((prevChromEnd, len(chrom.seq)))
    chromInsertions.append("")
    
    # Prior to masking the genome, create the list of insertions.
    chromSeq = chrom.seq.tostring()

    # now maskOut should be non-overlapping coordinates.
    modContig = ''.join([ chromSeq[chromIntervals[i][0]:chromIntervals[i][1]] + chromInsertions[i] for i in range(0,len(chromIntervals))])
    if (args.patchedGenome is not None):
        SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(modContig), id=chrom.id, name="", description=""), modGenome, "fasta")

    if (args.mode == 'patch'):
        # Produce a chain for this contig.
        #  score -- chain score
        #  tName -- chromosome (reference sequence)
        #  tSize -- chromosome size (reference sequence)
        #  tStrand -- strand (reference sequence)
        #  tStart -- alignment start position (reference sequence)
        #  tEnd -- alignment end position (reference sequence)
        #  qName -- chromosome (query sequence)
        #  qSize -- chromosome size (query sequence)
        #  qStrand -- strand (query sequence)
        #  qStart -- alignment start position (query sequence)
        #  qEnd -- alignment end position (query sequence)
        #  id -- chain ID
        chainFile.write("chain {} {} {} {} {} {} {} {} {} {} {} {}\n".format(1000, chrom.id, len(modContig), '+', 0, len(modContig), chrom.id, len(chrom), "+", 0, len(chrom), chainId))

        # Now produce description of alignment
        prevEnd = 0
        curPos = 0
        for i in range(0,len(patchesByChrom[name])):
            if (prevEnd > patchesByChrom[name][i][0]):
                pdb.set_trace()
                print "negative value"
            size = patchesByChrom[name][i][0] - prevEnd
            curPos += size
            if (args.liftedPatch is not None):
                liftBed.write("{}\t{}\t{}\t{}\n".format(name, curPos, curPos + len(patchesByChrom[name][i][2]), patchesByChrom[name][i][2]))
                
            tGap = len(patchesByChrom[name][i][2])
            curPos += tGap
            qGap = patchesByChrom[name][i][1] - patchesByChrom[name][i][0]
            chainFile.write("{} {} {}\n".format(size, tGap, qGap))
            prevEnd = patchesByChrom[name][i][1]
        # Last line has no gap elements, followed by a blank line
        chainFile.write("{}\n\n".format(len(chrom.seq) - prevEnd))



if (chainFile is not None):
    chainFile.close()

if (args.patchedGenome is not None):
    modGenome.close()
                            
    
    
    
    

