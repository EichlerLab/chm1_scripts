#!/usr/bin/env python
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import argparse
import matplotlib.lines as mlines
import math
import numpy as np
import tempfile
import subprocess
import os
import Tools
import pdb
#import colorbrewer
import colorsys
import sys

class Match:
    def __init__(self, qname, qfragstart, qfragend, qstrand, qstart, qend, qlen, tname, tstrand, tstart, tend, tlen, pctident):
        self.qname = qname
        self.qfragstart = qfragstart
        self.qfragend = qfragend
        self.qstrand = qstrand
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen
        self.tname = tname
        self.tstrand = tstrand
        self.tstart = tstart
        self.tend = tend
        self.tlen = tlen
        self.pctident = pctident
        self.rgb = ()
        

def IndexSeq(seq, k):
    index = {}
    seq = seq.upper()
    for i in range(len(seq)- k + 1):
        kmer = seq[i:i+k]
        if (kmer not in index):
            index[kmer] = []
        index[kmer].append(i)
    return index

def DotDirectional(query, targetIndex, k, d):
    query = query.upper()
    hits = []
    for i in range(len(query) - k + 1):
        word = query[i:i+k]
        if (word in targetIndex):
            for pos in targetIndex[word]:
                hits.append((i, pos, d))
    return hits

def DotPlot(query, target, k, forwardonly=False):
    targetIndex = IndexSeq(target.upper(),k)
    forward = DotDirectional(query, targetIndex, k, 0)
    if (forwardonly == False):
        revTarget = Seq.Seq(target).reverse_complement().tostring().upper()
        revTargetIndex= IndexSeq(revTarget, k)
        reverse = DotDirectional(query, revTargetIndex, k, 1)
        return forward + reverse
    else:
        return forward


def CountRegion(matches, qStart, qEnd):
    return sum( [ i[0] >= qStart and i[0] <= qEnd for i in matches ])
    

def SmartRegion(query, target, matches, k):
    np.set_printoptions(threshold=np.nan)
    counts = np.zeros(len(query))
    minTarget = (None,None)
    maxTarget = (None,None)
    if (len(matches) == 0):
        return (matches, (0, len(target)))
    
    if (len(query) < len(target)):
        return (matches, (0, len(query)))
    else:
        binSize = max(len(query)/200,20)
        nbins = (len(query) - len(target))/20 + 1
        counts = np.zeros(nbins)
        minCount = 0
        minRegion = 0

        for i in range(0,nbins):
            counts[i] = CountRegion(matches, i*binSize, i*binSize + len(target))
            if (counts[i] > minCount):
                minCount = counts[i]
                minRegion = (i*binSize, i*binSize + len(target))
        newMatches = []
        
        for m in matches:
            if (m[0] >= minRegion[0] and m[0] <= minRegion[1]):
                mp = (m[0] - minRegion[0], m[1], m[2])
                newMatches.append(mp)
        print "min region: " + str(minRegion)
        print "nu matches: " + str(len(matches))
        return (newMatches, (minRegion[0], min(minRegion[1], len(query))))


def SetTickLists(raggedTicks, raggedTickOffset, xoffset=0):
    nticks = sum(len(t) for t in raggedTicks)
    tickLayoutPos = [0]*nticks
    tickLayoutText = [""] *nticks
    idx = 0
    for i in range(len(raggedTicks)):
        for j in range(len(raggedTicks[i])):
            tickLayoutPos[idx] = raggedTickOffset[i][j]
            tickLayoutText[idx] = str(int(raggedTicks[i][j] + xoffset))
            idx+=1
    return tickLayoutPos, tickLayoutText

def ParseMatchLine(ml):
    mv = ml.split()
# parse this format
#    0                                 1                2 3 4     5   6    7    8     9 10   11   12          
#    scf7180000000013_1000_2000/0_1000 scf7180000000013 0 0 -5000 100 1000 2000 11963 0 1000 1000 20968
    qm = queryCoordRe.match(mv[0])
    if (qm is not None):
        g = qm.groups()
        qname = g[0]
        qfragstart = int(g[1])
        qfragend   = int(g[2])
    else:
        qname = mv[0]
        qfragstart = 0
        qfragend = int(mv[12])

    tname = mv[1]
    qstrand = int(mv[2])
    tstrand = int(mv[3])
    pctident = float(mv[5])
    tstart  = int(mv[6])
    tend    = int(mv[7])
    tlen    = int(mv[8])
    qstart  = int(mv[9])
    qend    = int(mv[10])
    qlen    = int(mv[11])
    return Match(qname, qfragstart, qfragend, qstrand, qstart, qend, qlen, tname, tstrand, tstart, tend, tlen, pctident)

def SetPlotCoordinates(match, querycs, targetcs, queryi, targeti):
    qi = queryi[match.qname]
    ti = targeti[match.tname]

    queryRangeStart= querycs[qi]
    queryRangeEnd  = querycs[qi+1]
    targetRangeStart = targetcs[ti]
    targetRangeEnd = targetcs[ti+1]
    # Shift the query coordinates to be relative to the 
    qStart = match.qfragstart + match.qstart
    qEnd   = match.qfragstart + match.qend
    
    # now handle strand orientation
    tStart = match.tstart
    tEnd   = match.tend
    if (match.tstrand == 1):
        tStart = match.tlen - match.tend
        tEnd   = match.tlen - match.tstart
        temp = qStart
        qStart = qEnd
        qEnd   = temp
    tStart += targetRangeStart
    tEnd   += targetRangeStart
    qStart += queryRangeStart
    qEnd   += queryRangeStart
    match.plTStart = tStart
    match.plTEnd   = tEnd
    match.plQStart = qStart
    match.plQEnd   = qEnd


    
queryCoordRe = re.compile('(.*)_(\d+)_(\d+)/.*')

ap = argparse.ArgumentParser(description = "make a dot plot")

ap.add_argument("--query", help="Read the query sequence.  Specify a read using bas:bas_fil_name.bas.h5:hole_number, or a region using region:chr:000-100 file.bam (two arguments)", required=True, nargs="+")
ap.add_argument("--target", help="Read the target sequence.", required=True)
ap.add_argument("--tappend", help="Also dot plot this sequence.", default=None)
ap.add_argument("--matches", help="Read the matches (m1 format), or dot:k to do a simple dot plot", required=False, default="dot:13")
ap.add_argument("--savefig", help="Write figure to a file.", default=None)
ap.add_argument("--region", help="Specify a region in the format ref:start-end .",default=None)
ap.add_argument("--regionFile", help="Specify a file with a region region in the format ref:start-end .",default=None)
ap.add_argument("--queryStart", help="Only use part of the query sequence starting here.", type=int, default=None)
ap.add_argument("--queryEnd", help="Use the part of the query sequence ending here.", type=int, default=None)
ap.add_argument("--blpath", help="Path to alternate blasr.", default="")
ap.add_argument("--nodel", help="Do not delete temporary files.", action='store_true', default=False)
ap.add_argument("--nolegend", help="Hide the legend.", action='store_true', default=False)
ap.add_argument("--xlabel", help="Override x label.", default=None)
ap.add_argument("--ylabel", help="Override y label.", default=None)
ap.add_argument("--xstart", help="Offset xtick label by this amount.", default=1, type=int)
ap.add_argument("--ystart", help="Offset ytick label by this amount.", default=1, type=int)
ap.add_argument("--title", help="use this as a title.", default=None)
ap.add_argument("--slop", help="Widen the region.", default=0, type=int)
ap.add_argument("--maxq", help="Maximum number of query sequences to align.  Randomly subsample if there are more queries than this.", default=None, type=int)
ap.add_argument("--forwardonly", help="Do not do reverse complement dotplot", default=False, action='store_true')
ap.add_argument("--smartadjust", help="Adjust the range on the query according to what likely aligned to the reference, trimming the ends.", action='store_true', default=False)
ap.add_argument("--thin", help="Draw thin lines", action='store_true', default=False)

args = ap.parse_args()


tmpFiles = []
if (args.ylabel is not None):
    ylabel = args.ylabel
else:
    ylabel = ""
xlabel = ""
# set up query information
queries = []
queryNames     = []

if (len(args.query) == 2 and args.query[0].find("region:") == 0):
    tmpQueryName = tempfile.mktemp(suffix=".query.fasta", dir=".")
    tmpFiles.append(tmpQueryName)
    command = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/RegionToFasta.py {} {} --out {}".format(args.query[1], args.query[0][7:], tmpQueryName)
    subprocess.call(command.split())
    queryIn = open(tmpQueryName)
    queryRecords = list(SeqIO.parse(queryIn, "fasta"))
    lengths = [len(rec.seq) for rec in queryRecords]

    for rec in queryRecords:
        queries.append(rec)
        queryNames.append(rec.id)
    queryIn.close()
                    
else:    
    for queryName in args.query:
        if (queryName[0:4] == "bas:"):
            # The input is a read, not a fasta file, extract that for the run.
        
            readValues = queryName.split(':')
            basName    = readValues[1]
            holeNumber = int(readValues[2])
            tmpQueryName = tempfile.mktemp(suffix=".query.fasta", dir=".")
            tmpFiles.append(tmpQueryName)
            command = "/net/eichler/vol5/home/mchaisso/software/blasr_1/cpp/pbihdfutils/bin/pls2fasta {} {} -holeNumber {} -trimByRegion".format(basName, tmpQueryName, holeNumber)
            args.query = tmpQueryName
            subprocess.call(command.split())
            ylabel = os.path.basename(basName) + " " + str(holeNumber)
            queryFileName = tmpQueryName
        else:
            queryFileName = queryName
            
        queryIn = open(queryFileName)
        queryRecords = list(SeqIO.parse(queryIn, "fasta"))
        lengths = [len(rec.seq) for rec in queryRecords]
        for rec in queryRecords:
            queries.append(rec)
            queryNames.append(rec.id)
        queryIn.close()
        
# set up taget
if (args.regionFile is not None):
    rf = open(args.regionFile)
    args.region = rf.readline().strip()
    
if (args.region is not None):
    faiFileName = args.target + ".fai"
    if (os.path.exists(faiFileName) == False):
        print "Error, when using a region there must be an associated .fai file."
        sys.exit(1)
    targetfd = file(args.target, 'r')

    fai = Tools.ReadFAIFile(faiFileName)
    region = Tools.ParseRegionStr(args.region)
    region2 = (region[0], max(0, region[1] - args.slop), min(fai[region[0]][0], region[2] + args.slop))
    print "REGION"
    print str(region2)
    seq = Tools.ExtractSeq(region2, targetfd, fai)
    tmpTarget = tempfile.NamedTemporaryFile(suffix=".target.fasta", dir=".", delete=False)
#    tmpTarget = open(tmpTargetName, 'w')
    tmpFiles.append(tmpTarget.name)
    record = SeqRecord.SeqRecord(Seq.Seq(seq), name="",id=args.region,description="")
    SeqIO.write(record, tmpTarget, "fasta")
    targetName = tmpTarget.name
    tmpTarget.close()
    args.target = targetName
    tmpFiles.append(targetName)    
    xlabel = args.region


#
# store information about targets
if (args.maxq is not None and args.maxq < len(queries)):
    r = np.random.permutation(len(queries))
    tmp = [queries[i] for i in r[0:args.maxq]]
    tmp2 = [queryNames[i] for i in r[0:args.maxq]]
    queries = tmp
    queryNames = tmp2
    
tf = open(args.target, 'r')
targets = list(SeqIO.parse(tf, "fasta"))

tf.close()

if (args.tappend is not None):
    tf = open(args.tappend)
    tseqs = list(SeqIO.parse(tf, "fasta"))
    targets.extend(tseqs)
    tf.close()
    xlabel += "    "
    xlabel += ' '.join([t.id for t in tseqs])
    
for q in range(0,len(queries)):
    queryStart = 0
    queryEnd = len(queries[q])
        

    if (args.queryEnd is not None):
        queryEnd = min(args.queryEnd, queryEnd)
    if (args.queryStart is not None):
        queryStart = min(args.queryStart, queryEnd)
    queries[q] = queries[q][queryStart:queryEnd]
    

queryDict = {q.name: len(q.seq) for q in queries }
targetDict = {t.name: len(t.seq) for t in targets }
queryNames = [q.name for q in queries]
targetNames = [t.name for t in targets]
numTargets = len(targetDict)
queryLengths = [len(q.seq) for q in queries]
targetLengths = [len(t.seq) for t in targets]

targetCumSum = [0]*(len(targets)+1)
for i in range(len(targets)):
    targetCumSum[i+1] = targetCumSum[i] + len(targets[i].seq)
queryCumSum = [0]*(len(queries) + 1)
for i in range(len(queries)):
    queryCumSum[i+1] = queryCumSum[i] + len(queries[i].seq)

totalQuery = sum([len(q.seq) for q in queries])
totalTarget = sum([len(t.seq) for t in targets])
    
if (ylabel == ""):
    ylabel = queries[0].id
if (xlabel == ""):
    xlabel = "  ".join([t.id for t in targets])

    
tFrac = [float(targetCumSum[i])/targetCumSum[-1] for i in range(len(targetCumSum)-1)]
qFrac = [float(queryCumSum[i])/queryCumSum[-1] for i in range(len(queryCumSum)-1)]




queryi = { queryNames[i] : i for i in range(len(queryNames))}
targeti = {targetNames[i] : i for i in range(len(targetNames))}
querycsd = { queryNames[i] : queryCumSum[i] for i in range(len(queryNames))}
targetcsd = { targetNames[i] : targetCumSum[i] for i in range(len(targetNames))}
matches = []

HSV_tuples = [(x*1.0/len(queries), 0.5, 0.5) for x in range(len(queries))]
RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

if (args.matches == None):
    # it is necessary to generate the matches, using blasr or a dot plot.
    targetName = args.target
    matchFileName = tempfile.mktemp(suffix=".m4", dir=".")
    matchFile = open(matchFileName, 'w')
    tmpFiles.append(matchFileName)
    args.matches = matchFileName
    blasrExe = args.blpath + "blasr"
    
    command = "{} {} {} -m 1 -out {} -bestn 100 -preserveReadTitle ".format(blasrExe, args.query, args.target, args.matches)
    print command
    subprocess.call(command.split())
    useCMap = True
    
elif (args.matches[0:4] == "dot:"):
    vals = args.matches.split(":")
    k = int(vals[1])
    useCMap = False


    querySeqs = Seq.Seq("")
    targetSeqs = Seq.Seq("")
    allMatches = []
    qFragStart = 0
    qFragEnd   = 0
    queryIndex = 0
    for query in queries:
        qFragEnd = qFragEnd + len(query) 
        targetFile = open(args.target, 'r')
        
        for target in targets:
            dotPlotMatches = DotPlot(query.seq.tostring(), target.seq.tostring(), k, args.forwardonly)
            if (args.smartadjust):
                (fixedDotPlotMatches, queryCoordinates) = SmartRegion(query.seq.tostring(), target.seq.tostring(), dotPlotMatches, k)

                args.ystart = queryCoordinates[0]
                query = query[queryCoordinates[0]:queryCoordinates[1]]

                dotPlotMatches = fixedDotPlotMatches
                queryCumSum = [0,len(query)]

            for m in dotPlotMatches:
                dpm = Match(query.id, 0, 0, 0, m[0], m[0]+k, len(query), target.id, m[2], m[1], m[1]+k, len(target), 100)
                dpm.rgb =RGB_tuples[queryIndex]
                SetPlotCoordinates(dpm, queryCumSum, targetCumSum, queryi, targeti)
                matches.append(dpm)
        queryIndex +=1
        targetFile.close()
        qFragStart = qFragEnd
        querySeqs = querySeqs + query

    args.matches = None
    
if (args.matches is not None):    
    mf = open(args.matches, 'r')
    for matchLine in mf:
        m = ParseMatchLine(matchLine)
        if (m is not None):
            SetPlotCoordinates(m, queryCumSum, targetCumSum, queryi, targeti)
            matches.append(m)

#
####################################################
# Set tick lists
####################################################
        
tStep = max(1,math.pow(10,math.floor(math.log10(targetCumSum[-1]))-1))

if (targetCumSum[-1]/tStep > 30):
    tStep = max(1,math.pow(10,math.floor(math.log10(targetCumSum[-1]))))
if (targetCumSum[-1]/tStep > 10):
    print "resetting ticks from " + str(tStep)
    tStep = int(np.floor(targetCumSum[-1]/10))
    print "to " + str(tStep)

qStep = max(1,math.pow(10,math.floor(math.log10(queryCumSum[-1]))-1))
if (queryCumSum[-1]/qStep > 10):
    qStep = int(np.floor(queryCumSum[-1]/10))

tTicks = [ np.arange(0, t, step=tStep) for t in targetLengths ]
qTicks = [ np.arange(0, q, step=qStep) for q in queryLengths]


tTickPositions = [ tTicks[i] + targetCumSum[i] for i in range(len(targetLengths))]
qTickPositions = [ qTicks[i] + queryCumSum[i] for i in range(len(queryLengths))]

ntTicks = sum(len(t) for t in tTicks)
nqTicks = sum(len(q) for q in qTicks)

            
tLayoutPos, tLayoutText = SetTickLists(tTicks, tTickPositions, args.xstart)
qLayoutPos, qLayoutText = SetTickLists(qTicks, qTickPositions, args.ystart)

lineWidth=2.0
if (args.thin):
    lineWidth=1.0
                
cmap = plt.get_cmap('jet')
gsi = 0
if (useCMap):
    rects = [mlines.Line2D([m.plTStart, m.plTEnd], [m.plQStart, m.plQEnd], color=cmap(1 + math.log(m.pctident/100.0)), lw=lineWidth) for m in matches]
    addLegend = False
else:
    rects = [mlines.Line2D([m.plTStart, m.plTEnd], [m.plQStart, m.plQEnd], color="black", lw=lineWidth) for m in matches]
    addLegend = True
    artists = [plt.Rectangle((0,0),1,1,color=RGB_tuples[i]) for i in range(len(queries))]
    

queryDim = 8.0
targetDim = queryDim * (max(1, 2*(float(totalTarget)/totalQuery))) + 2
#plt.figure(figsize=(targetDim, queryDim))
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8






plt.figure(figsize=(10,10), dpi=300)

if (len(queryNames) > 1):
    ax = plt.axes([0.07,0.15,0.40,0.80])
else:
    ax = plt.axes([0.2, 0.2, 0.70, 0.70])

for r in rects:
    ax.add_line(r)
    
ax.set_xticks(tLayoutPos)
ax.set_xticklabels(tLayoutText, rotation=90, size='x-large')
ax.set_yticks(qLayoutPos)
ax.set_yticklabels(qLayoutText, size='x-large')

#
# Attempt to get a good offset from the axes for where to write axis labels.
#
xSteps = np.floor(np.log10(args.xstart))
ySteps = np.floor(np.log10(max(queryLengths)))
print "steps: " + str(xSteps) + " " + str(ySteps)

ax.xaxis.set_label_coords(0.5, -0.15 + -0.01*xSteps)
ax.yaxis.set_label_coords(-0.15 + -0.01*ySteps, 0.5)

plt.xlabel(xlabel, fontsize=18)
if (ylabel != ""):
    plt.ylabel(ylabel, fontsize=18)

plt.subplots_adjust(bottom=0.3)

plt.xlim(xmin=0, xmax=targetCumSum[-1]*1.01)
plt.ylim(ymin=0, ymax=queryCumSum[-1]*1.01)

for t in range(len(targetCumSum)-1):
    plt.axvline(targetCumSum[t], linestyle='--', color='black')
for q in range(len(queryCumSum)-1):
    plt.axhline(queryCumSum[q], linestyle='--', color='black')
    
if (addLegend and args.nolegend == False):
    artists.reverse()
    revNames = queryNames
    revNames.reverse()
    plt.legend(artists, revNames, fontsize='small', shadow=True, bbox_to_anchor=(0,0,1,.95),  bbox_transform=plt.gcf().transFigure)


if (args.title is not None):
    plt.title(args.title)
        
if (args.savefig is not None):
    plt.savefig(args.savefig)
    cwd = os.getcwd()
    print "output is in " + cwd + "/" + args.savefig
else:
    plt.ioff()
    plt.show(block=True)
    
    
for i in range(len(targets)):
    plt.figtext(tFrac[i], 0.0, targets[0].name)




def get_open_fds():
    '''
    return the number of open file descriptors for current process

    .. warning: will only work on UNIX-like os-es.
    '''
    import subprocess
    import os

    pid = os.getpid()
    procs = subprocess.check_output( 
        [ "lsof", '-w', '-Ff', "-p", str( pid ) ] )

    procs =  filter(  lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(), procs.split( '\n' ) )
    fds = [int(a[1:]) for a in procs]
    return fds



if (args.nodel is False):
    openFds = get_open_fds()
    if (len(openFds) > 0):
        os.closerange(openFds[0], openFds[-1])
    for tmpFile in tmpFiles:
        command = "/bin/rm -f " + tmpFile
        try:
            print "command: " + command
            subprocesss.Popen(command.split())
            subprocess.wait()
        except:
            print "Cannot delete " + tmpFile
