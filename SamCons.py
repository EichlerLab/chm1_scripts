#!/usr/bin/env python


import argparse
import pysam
import sys
import Tools
import pdb

DELTA = 50

def Close(a, b):
    return (abs(a-b) < DELTA)
        
#    ---*==========*-----
#           --*===*----------
fromTo = 0
contains = 1
partial = 2

def Order(aln, boundaries):
    tb = boundaries[aln.tName]
    qb = boundaries[aln.title]
    min
    if (Close(tb[0], aln.tStart) and not Close(qb[0], aln.qStart)):
        return (aln.title, aln.tName, fromTo, aln.strand, aln.tStart)
    if (Close(tb[1], aln.tEnd) and not Close(qb[1], aln.qEnd)):
        return (aln.tName, aln.title, fromTo, aln.strand, aln.qStart)
    if (not Close(tb[0], aln.tStart) and not Close(tb[1], aln.tEnd) and
        Close(qb[0], aln.qStart) and Close(qb[1], aln.qEnd)):
        return (aln.tName, aln.title, contains, aln.strand, aln.tStart)
    if (Close(tb[0], aln.tStart) and Close(tb[1], aln.tEnd) and
        not Close(qb[0], aln.qStart) and not Close(qb[1], aln.qEnd)):
        return (aln.tName, aln.title, contains, aln.strand, aln.qStart)
    else:
        # partial alignment
        if (aln.tStart - tb[0] > aln.qStart - qb[0]):
            return (aln.tName, aln.title, partial, aln.tStart)
        else:
            return (aln.title, aln.tName, partial, aln.qStart)

class Vertex:
    def __init__(self, pos, read):
        self.pos    = [pos]
        self.read   = [read]
        self.outRed = []
        self.inRed  = []
        self.black  = [] # black edges are multi-directional
        self.index  = 0
        self.visited = 0

class Read:
    def __init__(self, name, vertexList):
        self.name = name
        self.vertexList = vertexList
        
def AssignIndices(reads):
    index = 0
    for read in reads.values():
        for vertex in read:
            vertex.index = index
            index += 1
            
def InitializeReadVertexList(length, seqIndex, k):
    vlist = [Vertex(i,seqIndex) for i in range(0,length-k+1)]
    for i in range(1,length-k):
        vlist[i-1].outRed.append(vlist[i])
        vlist[i-1].inRed.append(vlist[i])
    return vlist



ap = argparse.ArgumentParser(description="Build simple consensus from pairwise alignments.")
ap.add_argument("sam", help="Input sam file")
ap.add_argument("--k", help="Block size.", default=7, type=int)
ap.add_argument("--kgraph", help="Write the graph out.", default=None)
ap.add_argument("--rgraph", help="Write the reads-graph out.", default=None)
ap.add_argument("--minCount", help="Only print vertices if they have at least this many out edges.", default=0, type=int)
ap.add_argument("--delta", help="Wiggle room for alignment endpoints.", dest='DELTA')
args = ap.parse_args()

def SetBoundaries(samFileName, readBoundaries):
    samFile = open(samFileName, 'r')
    line = samFile.readline()
    while (line[0] == '@'):
        line = samFile.readline()
    while (line != ""):
        aln = Tools.SAMEntry(line)
        line = samFile.readline()
        if (aln.tName == aln.title):
            continue
        if (aln.tName not in readBoundaries):
            readBoundaries[aln.tName] = [aln.tStart, aln.tEnd]
        else:
            readBoundaries[aln.tName][0] = min(readBoundaries[aln.tName][0], aln.tStart)
            readBoundaries[aln.tName][1] = max(readBoundaries[aln.tName][1], aln.tEnd)
        if (aln.title not in readBoundaries):
            readBoundaries[aln.title] = [aln.qStart, aln.qEnd]
        else:
            readBoundaries[aln.title][0] = min(readBoundaries[aln.title][0], aln.qStart)
            readBoundaries[aln.title][1] = max(readBoundaries[aln.title][1], aln.qEnd)

    
def StoreReadGraph(samFileName, readBoundaries, readGraph):
    samFile = open(samFileName, 'r')
    line = samFile.readline()
    while (line[0] == '@'):
        line = samFile.readline()
    while (line is not None and line != ""):
        aln = Tools.SAMEntry(line)
        if (aln.tName == aln.title):
            line = samFile.readline()
            continue
        order = Order(aln, readBoundaries)
        if (order is None):
            print "t: {},{}: {}-{}  q: {},{}: {}-{} strand {}".format(readBoundaries[aln.tName][0],readBoundaries[aln.tName][1],
                                                            aln.tStart, aln.tEnd,
                                                            readBoundaries[aln.title][0],readBoundaries[aln.title][1],
                                                            aln.qStart, aln.qEnd, aln.strand)
        else:
            if (order[0] not in readGraph):
                readGraph[order[0]] = {}
            if (order[1] not in readGraph[order[0]]):
                readGraph[order[0]][order[1]] = order[2]
        line = samFile.readline()
            
def MakeForward(pos, strand, length, k):
    if (strand == 0):
        return pos
    else:
        return length - pos - k


def PrintGraph(reads, readNames, graphName):
    graphOut = open(graphName, 'w')
    # first print read edges.
#    graphOut.write("graph qvgraph {\n")
    for read in reads.values():
        for i in range(len(read)-1):
            if (len(read[i].black) >= args.minCount and
                len(read[i+1].black) >= args.minCount):
                graphOut.write("  {} red {} \n".format(read[i].index, read[i+1].index))
    # now output black edges
    for name,read in reads.iteritems():
        for i in range(len(read)):
            if (len(read[i].black) < args.minCount):
                continue
            for blackEdge in read[i].black:
                destRead = reads[readNames[blackEdge[0]]]
                if (name > destRead):
                    pass
                if (len(destRead[blackEdge[1]].black) >= args.minCount):
#                    graphOut.write("  {} -> {} [ color=\"black\" ]\n".format(read[i].index, destRead[blackEdge[1]].index))
                    graphOut.write("  {} black {} \n".format(read[i].index, destRead[blackEdge[1]].index))
    graphOut.close()
#    graphOut.write("}\n")


    

samFile= open(args.sam, 'r')
line = samFile.readline()
reads = {}
readIndex = 0
readLengths = {}
readIndex = {}
index = 0
readNames = []

readBoundaries = {}
SetBoundaries(args.sam, readBoundaries)
for k,v in readBoundaries.iteritems():
    print k + " " + str(v)

readGraph = {}
print "storing read graph."
StoreReadGraph(args.sam, readBoundaries, readGraph)
if (args.rgraph is not None):
    graphOut = open(args.rgraph, 'w')
    for src in readGraph.keys():
        for dest in readGraph[src].keys():
            if (readGraph[src][dest] != contains):
                graphOut.write(str(src) + " " + str(readGraph[src][dest]) + " " + str(dest) + "\n")
    graphOut.close()


while (line[0] == '@'):
    if (line[0:3] == "@SQ"):
        vals = line[3:].split()
        name = Tools.GetKV("SN:", vals)
        length = int(Tools.GetKV("LN:", vals))
        reads[name] = InitializeReadVertexList(length, index, args.k)
        print "initialized " + name
        readLengths[name] = length
        readIndex[name] = index
        readNames.append(name)
        index +=1
    line = samFile.readline()
#
# Make a list that may be indexed by number and not name.
#
readList = [ Read(name, reads[name]) for name in readNames ]
prevQuery = ""
while (line != ""): 
    aln = Tools.SAMEntry(line)
    qPos = aln.qStart - 1
    tPos = aln.tStart - 1
    tList = reads[aln.tName]
    qList = reads[aln.title]
    tLength = len(tList)+args.k-1
    qLength = len(qList)+args.k-1
    tIndex = readIndex[aln.tName]
    qIndex = readIndex[aln.title]

    q = 0
    t = 0
    qForward = 0
    if (aln.tName == aln.title):
        line = samFile.readline()
        print "skipping  () " + aln.title + " " + aln.tName
        continue
    #
    # Apply some filtering for which alignments are considered.
    #
    # rule 1. target must be longer than the query
    if (len(tList) < len(qList)):
        
        continue
    # rule 2. only one alignment to the target per query
    if (aln.title != prevQuery):
        targetAlignCount = {}
    else:
        if (aln.tName not in targetAlignCount):
            targetAlignCount[aln.tName] = 1
        else:
            print "skipping (multi) " + aln.title + " " + aln.tName
            continue
    prevQuery = aln.title
    
    nMatches = 0
    for i in range(len(aln.ops)):
        if (aln.ops[i] == 'M'):
            # the alignment contains a sufficient match. 
            if (aln.lengths[i] >= args.k):
                for i1 in range(aln.lengths[i] - args.k + 1):
                    q = qPos + i1
                    t = tPos + i1
                    qForward = MakeForward(q, aln.strand, qLength, args.k)
                    if (t > len(tList)):
                        print "error! " + str(i) + " " + str(len(aln.ops))+ " " + str(t) + " " + str(len(tList))
                        sys.exit(1)
                    if (qForward >= len(qList)):
                        print "error! " + str(i) + " " + str(qForward) + " " + str(len(qList))
                        sys.exit(1)
                    tList[t].black.append((qIndex, qForward, aln.strand))
                    qList[qForward].black.append((tIndex, t, aln.strand))
                    nMatches += 1
            qPos += aln.lengths[i]
            tPos += aln.lengths[i]
        elif (aln.ops[i] == 'D'):
            tPos += aln.lengths[i]
        elif (aln.ops[i] == 'I'):
            qPos += aln.lengths[i]
    line = samFile.readline()

def AddVertexToConflictSet(vertex, conflictSet):
    i = vertex.read[0]
    if (i not in conflictSet):
        conflictSet[i] = []
    conflictSet[i].append(vertex.pos[0])


def SearchForConflict(readList, vertex, conflictSet, level=1):
    AddVertexToConflictSet(vertex, conflictSet)
    vertex.visited = 1
#    print str(vertex.read)
#    print str(vertex.pos)
#    print "level: " + str(level)

    
    for outEdge in vertex.black:
#        print str(level) + " visiting " + str(outEdge) 
        if (readList[outEdge[0]].vertexList[outEdge[1]].visited == 0):
            readList[outEdge[0]].vertexList[outEdge[1]].visited = 1
            SearchForConflict(readList, readList[outEdge[0]].vertexList[outEdge[1]], conflictSet, level+1)
    print "done visiting all out edges for " + str(vertex)

def IsVertexConflicted(readList, conflictSet):
    for k,v in conflictSet.iteritems():
        if (len(v) > 1):
            print k
            print v
            print readList[k].name + " " + ', '.join([str(a) for a in v])
            return True
    return False
for vlist in reads.values():
    for i in range(len(vlist)):
        if (len(vlist[i].black) > 0):
            conflictSet = {}
            SearchForConflict(readList, vlist[i], conflictSet)
            if (IsVertexConflicted(readList, conflictSet)):
                print "found conflict"
                print str(conflictSet)
            print "Searched for conflicts for " + str(i)
#for name in reads.keys():
#    i = 0;
#
#    prevPos = 0
#    vlist = reads[name]
#    nVertex = len(vlist)
#    steps = []
#    while (i < nVertex):
#        while (i < nVertex and len(vlist[i].black) == 0):
#            i += 1
#        if (i < nVertex):
#            if (len(vlist[i].black) == 1):
#                other = vlist[i].black[0][0]
#            else:
#                other = "*"
#            steps.append((i - prevPos, len(vlist[i].black), other))
#        prevPos = i
#        i+=1
#    steps.append((i - prevPos, 0, "*"))
#    print str(steps)

AssignIndices(reads)

if (args.kgraph is not None):
    PrintGraph(reads, readNames, args.kgraph)
