#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
import tempfile
import argparse
import subprocess

ap = argparse.ArgumentParser(description="Find inverted duplications in two sequences.")

ap.add_argument("query", help="Query sequence.")
ap.add_argument("target", help="Target sequence.")
ap.add_argument("--region", help="Use this region of the target.", default=None)
ap.add_argument("--slop", help="Slop", default=0, type=int)
ap.add_argument("--k", help="word size", type=int, default=9)
ap.add_argument("--min", help="Min anchor size", type=int, default=5)
ap.add_argument("--gap", help="Max separating gap", type=int, default=20)

args = ap.parse_args()


import Tools
(chrom,start,end) = Tools.ParseRegionStr(args.region)

# Will need to make a file for the reference
start = max(0, start - args.slop)
region = "{}:{}-{}".format(chrom,start, end+args.slop)

tempRef = None
if (args.region is not None):
    tempRef = tempfile.NamedTemporaryFile(suffix=".fasta")
    print tempRef.name
    command = "samtools faidx {} {}".format(args.target, region)
    subprocess.call(command.split(), stdout=tempRef)
    args.target = tempRef.name
    
querySeq = SeqIO.read(open(args.query), "fasta")
targetSeq = SeqIO.read(open(args.target), "fasta")


def StoreUniqueWords(seq, k):
    allWords = {}
    for i in range(0,len(seq)-k+1):
        sseq = seq[i:i+k]
        if (sseq in allWords):
            allWords[sseq] = False
        else:
            allWords[sseq] = i
    uniqueWords = {}
    for k in allWords.keys():
        if (allWords[k] != False):
            uniqueWords[k] = allWords[k]
    return uniqueWords



uniqueTarget = StoreUniqueWords(str(targetSeq.seq), args.k)


def Matches(seq, words, k):
    rcMatches = []
    for i in range(0,len(seq)-k+1):
        sstr = seq[i:i+k]
        if (sstr in words):
            rcMatches.append((i,words[sstr]) )
    return rcMatches
            


ftm = Matches(str(querySeq.seq), uniqueTarget, args.k)
#print ftm[0:10]
#print ftm[-10:-1]



def MergeMatches(m):
    matches = []
    i = 0
    curMatch = [m[i][0], m[i][1], m[i][0], m[i][1]]
    nMatches = 1
    while (i + 1 < len(m)):
        da = m[i+1][0] - m[i][0]
        db = m[i+1][1] - m[i][1]
        if (da - db < 4 and m[i+1][0] - m[i][0] < args.gap ):
            curMatch[2] = m[i+1][0]
            curMatch[3] = m[i+1][1]
            nMatches += 1
        else:
            if (nMatches >= args.min):
                matches.append((curMatch, nMatches))
            curMatch = [m[i+1][0], m[i+1][1], m[i+1][0], m[i+1][1]]
            nMatches = 1
        i += 1
    curMatch[2] = m[i][0]
    curMatch[3] = m[i][1]
    if (nMatches >= args.min):
        matches.append((curMatch, nMatches) )
    return matches
        
orthMatches =  MergeMatches(ftm)

queryStart = orthMatches[0][0][0]
queryEnd   = orthMatches[-1][0][2]
targetStart = orthMatches[0][0][1]
targetEnd   = orthMatches[-1][0][3]

print "{} {} {} {}".format(queryStart, queryEnd, targetStart, targetEnd)

uniqueRevQuery = StoreUniqueWords(str(querySeq.seq.reverse_complement()), args.k)
uniqueRevTarget = StoreUniqueWords(str(targetSeq.seq.reverse_complement()), args.k)
queryRevComp = str(querySeq.reverse_complement())

qm = Matches(str(querySeq.seq), uniqueRevTarget, args.k)
tm = Matches(str(targetSeq.seq), uniqueRevQuery, args.k)

print "{} {}".format(len(uniqueRevQuery), len(uniqueRevTarget))
#print qm    
mqm = MergeMatches(qm)
print "against rev target:"
for i in range(0,len(mqm)):
    qs = mqm[i][0][0]
    qe = mqm[i][0][2]
    ts = len(targetSeq.seq) - mqm[i][0][3]
    te = len(targetSeq.seq) - mqm[i][0][1]
    if (qs >= queryStart and qe <= queryEnd and ts >= targetStart and te <= targetEnd):
        print "{} {} {} {} {}".format(qs, qe, ts, te, mqm[i][1])

print "against rev query:"
mtm = MergeMatches(tm)

for i in range(0,len(mtm)):
    qs = len(querySeq.seq) - mtm[i][0][2]
    qe = len(querySeq.seq) - mtm[i][0][1]
    ts = mtm[i][0][1]
    te = mtm[i][0][3]
#    if (qs >= queryStart and qe <= queryEnd and ts >= targetStart and te <= targetEnd):
    print "{} {} {} {} {}".format(qs, qe, ts, te, mqm[i][1])
