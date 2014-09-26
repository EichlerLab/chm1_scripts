#!/usr/bin/env python

import pysam
import argparse
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import pdb
import bisect

ap = argparse.ArgumentParser(description="Cluster breakpoints.")
ap.add_argument("table", help="Input tabular file")
ap.add_argument("left", help="Left clusters.", nargs=2)
ap.add_argument("right", help="Right clusters.", nargs=2)
ap.add_argument("out", help="Output file, - implies stdout")
ap.add_argument("--delta", help="Allowable gap between endponts.", type=int, default=1000)
ap.add_argument("--minSupport", help="Min overlapping clusters", type=int, default=2)
args = ap.parse_args()

tabFile = open(args.table)
#
#   *-----*
#     *---*  right-overlap
#  *---*     left-overlap, and close enough

def DefineClusters(table, span, col, minSupport = 2):
    i = 0
    j = 1
    clusters = {}
    while ( i < len(table) and j < len(table) ):
        j = i + 1
        if (j == len(table)):
            break
        maxSpan = table[j][2] - table[j][1]
        maxSpanIndex = i
        minLeft = table[i][1]
        maxRight = table[i][2]
        while (j < len(table) and
               table[j][0] == table[j-1][0] and
               table[j][col] < table[j-1][col] + span):
            if (minLeft > table[j][1]):
                minLeft = table[j][1]
            if (maxRight < table[j][1]):
                maxRight = table[j][2]
            curSpan = maxRight - minLeft
            if (curSpan > maxSpan):
                maxSpan = curSpan
                maxSpanIndex = j

            j += 1
                
        if (table[i][0] not in clusters):
            clusters[table[i][0]] = []
        if (j - i >= minSupport):
            clusters[table[i][0]].append((minLeft, maxRight, j - i))
        i = j
    
    return clusters


def ReadClusters(clusterFileName, side, minSupport=2, span=100):
    cf = open(clusterFileName, 'r')
    clusters = {}
    nAdded   = 0
    table = []
    for line in cf:
        vals    = line.split()
        chrom   = vals[0]
        start   = int(vals[1])
        end     = int(vals[2])
        table.append((chrom,start,end))

    return DefineClusters(table, span, side, minSupport)
    
def FindCluster(chrom, start, clusters):
    s = (start, start)
    if (chrom not in clusters):
        return None
    sb = bisect.bisect_left(clusters[chrom], s)

    while (sb < len(clusters[chrom]) and sb > 0 and clusters[chrom][sb][0] > start):
        sb -= 1
    if (sb >= len(clusters[chrom])):
        return None
    if (clusters[chrom][sb][0] <= start and clusters[chrom][sb][1] >= start):
        return (chrom, sb)
    return None
    

class Interval:
    def __init__(self):
        self.intv = []

    def Check(self, lp, rp):
        for i in range(len(self.intv)):
            if (abs(lp - self.intv[i][0]) < args.delta or
                abs(rp - self.intv[i][1]) < args.delta):
                return True
        return False

    def Add(self, lp, rp):
        self.intv.append((lp,rp))
#        self.intv.sort()

    def Size(self):
        return len(self.intv)

    def Coordinates(self):
        return (np.min(self.intv) , np.max(self.intv))
        
class Cluster:
    def __init__(self, leftChr, leftStrand, lstart, lend, rightChr, rightStrand, rstart, rend):
        self.leftStrand = leftStrand
        self.rightStrand = rightStrand
        self.leftChr = leftChr
        self.rightChr = rightChr
        self.leftBp = Interval()
        self.rightBp = Interval()
        if (lend< rend):
            self.leftBp.Add(lstart, lend)
            self.rightBp.Add(rstart, rend)
        else:
            self.leftBp.Add(rstart, rend)
            self.rightBp.Add(lstart, lend)
        
    def CoordinatesOverlap(self, left, right, side):
        return getattr(self, side).Check(left, right)

    def IntervalsOverlap(self, lchr, lstr, lsta, lend, rchr, rstr, rsta, rend):
        if ((lchr == self.leftChr and lstr == self.leftStrand and self.CoordinatesOverlap(lsta, lend, 'leftBp')) and
            (rchr == self.rightChr and rstr == self.rightStrand and self.CoordinatesOverlap(rsta, rend, 'rightBp'))):
            return 1
        elif ((lchr == self.rightChr and lstr == self.rightStrand and self.CoordinatesOverlap(lsta, lend, 'rightBp')) and
              (rchr == self.leftChr and  rstr == self.leftStrand  and self.CoordinatesOverlap(rsta, rend, 'leftBp'))):
            return 2
        else:
            return 0

    def Add(self, left, right, side):
        getattr(self, side).Add(left, right)
        

nLines = 0
if (args.out ==  "-"):
    outFile = sys.stdout
else:
    outFile = open(args.out, 'w')

allClusters = {}
leftClusters = {}
rightClusters = {}

llClusters = ReadClusters(args.left[0], 1, args.minSupport, args.delta)
lrClusters = ReadClusters(args.left[1], 2, args.minSupport, args.delta)
rlClusters = ReadClusters(args.right[0], 1, args.minSupport, args.delta)
rrClusters = ReadClusters(args.right[1], 2, args.minSupport, args.delta)

clusterPairCounts = {}
leftClusters  = [llClusters, lrClusters]
rightClusters = [rlClusters, rrClusters]

for line in tabFile:
    vals = line.split()
    lscore = int(vals[1])
    lstrand = int(vals[2])
    lchr = vals[3]
    lstart = int(vals[4])
    lend   = int(vals[5])
    rscore = int(vals[6])
    rstrand = int(vals[7])
    rchr = vals[8]
    rstart = int(vals[9])
    rend   = int(vals[10])


    clLeftStart  = FindCluster(lchr, lstart, llClusters)
    clLeftEnd    = FindCluster(lchr, lend, lrClusters)
    clRightStart = FindCluster(rchr, rstart, rlClusters)
    clRightEnd   = FindCluster(rchr, rend, rrClusters)
#    print "coords"
#    print str((lstart, lend, rstart, rend))
#    print str((clLeftStart, clLeftEnd, clRightStart, clRightEnd))
    if ( (clLeftStart != None or clLeftEnd != None) and
         (clRightStart != None or clRightEnd != None) ):
        if (clLeftStart != None):
            leftCluster = (clLeftStart[0], clLeftStart[1], 0)
        else:
            leftCluster = (clLeftEnd[0],clLeftEnd[1], 1)
            
        if (clRightStart != None):
            rightCluster = (clRightStart[0], clRightStart[1], 0)
        else:
            rightCluster = (clRightEnd[0], clRightEnd[1], 1)

        if (leftCluster not in clusterPairCounts):
            clusterPairCounts[leftCluster] = {}
        if (rightCluster not in clusterPairCounts[leftCluster]):
            clusterPairCounts[leftCluster][rightCluster] = 0
        clusterPairCounts[leftCluster][rightCluster] += 1

for cpLeft in clusterPairCounts.keys():
    for cpRight in clusterPairCounts[cpLeft].keys():
        outFile.write(cpLeft[0] + "\t" + str(leftClusters[cpLeft[2]][cpLeft[0]][cpLeft[1]][0]) + "\t" + str(leftClusters[cpLeft[2]][cpLeft[0]][cpLeft[1]][1]) + "\t" + cpRight[0] + "\t" + str(rightClusters[cpRight[2]][cpRight[0]][cpRight[1]][0]) + "\t" + str(rightClusters[cpRight[2]][cpRight[0]][cpRight[1]][1]) + "\t" + str(clusterPairCounts[cpLeft][cpRight]) + "\n" )

if (outFile != sys.stdout):
    outFile.close()
    
##except:
##    print "problem"
