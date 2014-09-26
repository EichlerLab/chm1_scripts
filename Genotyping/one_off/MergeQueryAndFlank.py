#!/usr/bin/env python

import sys

if (len(sys.argv) < 5):
    print "usage: Merge4QueryAndFlank.py file.query flank.query window flank.fasta.gc output.query"
    sys.exit(0)
    
localFile = open(sys.argv[1])
flankingFile = open(sys.argv[2])
window = int(sys.argv[3])
flankLengths = open(sys.argv[4])
outputFile = open(sys.argv[5], 'w')

flankLengthLines = [line.split() for line in flankLengths]
flankLengths = { f[0]: int(f[3]) for f in flankLengthLines}


local = [ line.split() for line in localFile]
flanking = [line.split() for line in flankingFile]

localQuery = {}
for i in range(0,len(local)):
    if (local[i][0] not in localQuery):
        localQuery[local[i][0]] = {}
    localQuery[local[i][0]][local[i][1]] = True


#
# now merge files
#


curRegion = ""
li = 0
fi = 0
numRemoved = 0
for fi in range(0,len(flanking)):
    if (flanking[fi][0] not in flankLengths):
        print "ERROR, missing region: " + flanking[fi][0]
        sys.exit(0)
        
    regionLength = flankLengths[flanking[fi][0]]
    regionPos = int(flanking[fi][2])
    if (regionPos < window or regionPos + window > regionLength):
        outputFile.write( "\t".join(flanking[fi]) + "\n")
    else:
        if (flanking[fi][1] in localQuery[flanking[fi][0]]):
            outputFile.write("\t".join(flanking[fi]) + "\n")
        else:
            numRemoved += 1

sys.stderr.write("Removed " + str(numRemoved) + "\n")            

outputFile.close()
