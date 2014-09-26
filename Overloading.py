#!/usr/bin/env python


import sys
import re
import numpy
import argparse


subreadCoordsRe = re.compile('.*\/(\d+)_(\d+)')
kvpRe = re.compile('(.*):.:(.*)')


clp = argparse.ArgumentParser(description='Check for overloading of PacBio SMRTCells.')
clp.add_argument("--sam", help="Alignment file.")
clp.add_argument("--verbose", help="Print verbose summaries of bad reads.", action="store_true", default=False)

argVals = clp.parse_args()

samFile = open( argVals.sam , "r" )

while(samFile):
    line = samFile.readline()
    if (len(line) == 0):
        sys.exit(1)
    if (line[0] != '@'):
        break

    
prevLines = [line]
prevLine = line
prevVals = [line.split()]

for line in samFile.readlines():
    vals = line.split()
    if (vals[0] != prevVals[0][0]):
        #
        # handle multiple possible hits
        #

        status = "OK"


        m = subreadCoordsRe.match(prevVals[0][0])

        g = m.groups()
        if (len(g) != 2):
            print "ERROR matching " + prevVals[0][0]
        start = int(g[0])
        end   = int(g[1])
            #               if (len(prevLines) > 1):
        subreadLength = end - start
#        print str(prevLines)

        if (m):
            #                print prevVals[0]
            #                 print str((start,end))
            matches = [0]*subreadLength
            index = 0
            for aln in prevLines:
                index += 1
                alnVals = aln.split()
                cigar = alnVals[5]
                flag  = int(alnVals[1])
                xs = start
                xe = end
                for kvp in alnVals[11:]:
                    m = kvpRe.match(kvp)

                    if (m is not None):
                        g = m.groups()
                    if (g[0] == "XS"):
                        xs = int(g[1])
                    if (g[0] == "XE"):
                        xe = int(g[1])
                
                cigarOps = re.sub('[0-9]', ' ', cigar).split()
                cigarLens = [int(i) for i in re.sub('[A-Z]', ' ', cigar).split()]

                if (flag & 16 != 0):
                    cigarOps.reverse()
                    cigarLens.reverse()
                pos = xs - 1 - start
                #                    print "starting at " + str(pos)
                for i in range(0,len(cigarOps)):
                    if (cigarOps[i] == 'I'):
                        pos += cigarLens[i]
                    if (cigarOps[i] == 'M'):
                        for mi in range(pos, pos + cigarLens[i]):
#                            print str(matches)
                            if (mi < len(matches) and matches[mi] == 0):
                                matches[mi] = index
                        pos += cigarLens[i]
                        

            pos1, = numpy.where(numpy.diff(matches) != 0)
            pos1 = numpy.concatenate(([0],pos1+1,[len(matches)]))

            matchesRle = [[matches[a], b-a] for (a,b) in zip(pos1[:-1],pos1[1:])]

            runs = []
            i = 0
            while (i < len(matchesRle)):
                if (matchesRle[i][1] < 5):
                    del matchesRle[i]
                else:
                    i += 1

            for m in matchesRle:
                l = len(runs)
                if (l == 0):
                    if (m[0] != 0):
                        runs.append(m)
                else:
                    if (runs[l-1][0] == m[0]):
                        runs[l-1][1] += m[1]
                    else:
                        if (m[0] != 0):
                            runs.append(m)

            alns = [False]*len(prevLines)
            nRuns = 0
            for r in runs:
                if (r[1]> 100):
                    alns[r[0]-1] = True
                    nRuns += 1


            if (nRuns > 1):
                nAlns = len(prevLines)
                sameChr = [[False]*nAlns]*nAlns
                alnDiff = [[0]*nAlns]*nAlns

                for i in range(0,nAlns-1):
                    for j in range(1,nAlns):
                        sameChr[i][j] = prevVals[i][2] == prevVals[j][2]
#                        print str(i) + " " + str(j) + " " + str(prevVals[i][3]) + " " + str(prevVals[j][3]) + " " + str(abs(int(prevVals[i][3]) - int(prevVals[j][3])))
                        alnDiff[i][j] = abs(int(prevVals[i][3]) - int(prevVals[j][3]))

                for i in range(0,nAlns-1):
                    for j in range(1,nAlns):
                        if (alns[i] == False or alns[j] == False):
                            continue
                        if (sameChr[i][j] and alnDiff[i][j] < end - start):
                            
                            status = "MissingAdapter"
                        elif (sameChr[i][j] == False or (alnDiff[i][j] > 100000)):
                            status = "Chimera"
#                            print str(runs)
#                            print alnDiff[i][j]
        if (status != "OK" and argVals.verbose == True):
            print str(runs)
        print prevVals[0][0] + " " + str(subreadLength) + " " + status
        
        prevLines = [line]
        prevVals = [line.split()]
    else:
        prevLines.append(line)
        prevVals.append(line.split())
    prevLine = line

