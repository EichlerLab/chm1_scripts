#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import os
import sys


regionFile = open(sys.argv[1])
dirName    = sys.argv[2]
outputName = sys.argv[3]

outputFile = open(outputName, 'w')

def LineToRegion(line):
    vals = line.split()
    if (vals[0].find(":") >= 0):
        return vals[0]
    else:
        return vals[0]+":"+vals[1]+"-"+vals[2]

regionIndex = 1    
for line in regionFile:
    region = LineToRegion(line)
#    fastaName= "rgn_"+str(regionIndex) + "/" + dirName + "/9-terminator/"+dirName + ".ctg.fasta"
    fastaName= "rgn_"+str(regionIndex) + "/reads.consensus.fasta"
    if (os.path.exists(fastaName)):
        print fastaName
        fastaFile = open(fastaName)
        contigIndex = 0
        for seq in SeqIO.parse(fastaFile, "fasta"):
            seq.id = region + "/" + str(contigIndex)
            SeqIO.write(seq, outputFile, "fasta")
    regionIndex += 1
outputFile.close()
            
    
