#!/usr/bin/env python

import sys
import Tools

if (len(sys.argv) != 2):
    print "usage ScreenForGaps.py sam"
    print " sam.gaps.bed must exist."
    exit(1)
    
samFileName = sys.argv[1]
bedFileName = samFileName + ".gaps.bed"

samFile  = open(samFileName)
bedFile = open(bedFileName)

flank = 4000

mapping = {}
for line in samFile:
    if (line[0] == '@'):
        continue
    vals = line.split()
    samEntry = Tools.SAMEntry(line)
    nMatch = 0
    for i in range(0,len(samEntry.ops)):
        if (samEntry.ops[i] == 'M'):
            nMatch += samEntry.lengths[i]
    mapping[vals[0]] = [vals[2], int(vals[3]), nMatch]
    

prevQuery = ""
gapBases = 0
queryAlnBases = 0
print "query\tqueryGapBases\tgapBases\tqueryMatch\tqueryLen"
for line in bedFile:
    vals = line.split()
    gapChrom = vals[0]
    gapStart = int(vals[1])
    query = vals[-1]

    if (query != prevQuery and prevQuery != ""):
        print prevQuery + "\t" + str(queryGapBases) + "\t" + str(gapBases) + "\t" + str(queryAlnBases) + "\t" + str(seqLen)
        gapBases = 0

    if (query not in mapping):
        continue

    queryAlnBases = mapping[query][2]
    
    if (mapping[query][0] != gapChrom):
        continue


    # parse the seq len from the seq name
    seqLen = int(query.split("_")[-1])
    mappingStart = mapping[query][1]
    queryGapBases = seqLen - 2*flank;
    
    relGapStart = gapStart - mappingStart

    gapLen = int(vals[4])
    if ((relGapStart >= flank and relGapStart < seqLen - flank) or
        (relGapStart <= flank and relGapStart + gapLen > flank)):
        gapBases += int(vals[4])


    prevQuery = query

        
        
print query + "\t" + str(queryGapBases) + "\t" + str(gapBases) + "\t" + str(queryAlnBases) + "\t" + str(seqLen)

        
    
    
