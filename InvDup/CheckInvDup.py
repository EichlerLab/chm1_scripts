#!/usr/bin/env python


import argparse
import tempfile
import sys
import subprocess
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import os
import Tools

ap = argparse.ArgumentParser(description="Look to see if an assembly contains an inverted duplication")

ap.add_argument("assembly", help="Assembled sequence.")
ap.add_argument("--region", help="Region that was assembled.", default=None)
ap.add_argument("--regionFile", help="File with region that was assembled.", default=None)
ap.add_argument("--slop", help="Window around the region to search", default=None)
ap.add_argument("--reference", help="Specify where to find the reference", default="/net/eichler/vol2/eee_shared/assemblies/hg19/ucsc.hg19.fasta")
ap.add_argument("--summary", help="Write a summary to this file.", default=None)
ap.add_argument("--keeptemp", help="Keep temporary files around after run.", default=True, action='store_false', dest='deletetemp')
ap.add_argument("--deletion", help="Check for deletion by swapping reference and assembly.", default=False, action='store_true')


args = ap.parse_args()
if (args.regionFile is not None):
    rf = open(args.regionFile)
    args.region = rf.readline().strip()
    rf.close()

(chrom,start,end) = Tools.ParseRegionStr(args.region)

if (args.slop is None):
    fileSize = os.path.getsize(args.assembly)
    args.slop = int(fileSize / 2)

summaryFile = None    
if (args.summary is not None):
    summaryFile = open(args.summary, 'a+')
    
# Will need to make a file for the reference
start = max(0, start - args.slop)
region = "{}:{}-{}".format(chrom,start, end+args.slop)

refSeq = tempfile.NamedTemporaryFile(delete=args.deletetemp, suffix=".fasta")
samFile = tempfile.NamedTemporaryFile(delete=args.deletetemp, suffix=".sam")
gapBedFile = tempfile.NamedTemporaryFile(delete=args.deletetemp, suffix=".bed")

command = "samtools faidx {} {}".format(args.reference, region)
subprocess.call(command.split(), stdout=refSeq)

targetFile = refSeq.name
queryFile  = args.assembly
if (args.deletion):
    targetFile = args.assembly
    queryFile = refSeq.name

command = "blasr {} {} -affineAlign -affineOpen 100 -affineExtend 0 -sam -out {} -bestn 1".format(queryFile, targetFile, samFile.name)
subprocess.call(command.split())

command = "samtools faidx {}".format(targetFile)
subprocess.call(command.split())

command = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py {} {} --tsd 20 --condense 20 --outFile {} ".format(targetFile, samFile.name, gapBedFile.name)
subprocess.call(command.split())

command = "/bin/rm {}.fai".format(targetFile)
subprocess.call(command.split())

gapBedFile = open(gapBedFile.name)


print "START"
print args.assembly
print region
for line in gapBedFile:
    vals = line.split()

    gapFile = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    seqName = vals[0]+"/"+vals[1]+"/"+vals[2]
    rec = SeqRecord.SeqRecord(Seq.Seq(vals[5]), id=seqName, name="", description="")
    SeqIO.write(rec, gapFile, "fasta")
    gapFile.close()
    gapAlnFile = tempfile.NamedTemporaryFile(delete=args.deletetemp,suffix=".m1")

    if (vals[3] == "insertion"):
        target = queryFile
    else:
        target = targetFile
        
    #
    # with the insertion, map to the assembly for the possible inverted duplication
    #

    remapCommand = "blasr {} {}  -bestn 2 -sdpTupleSize 11 -m 1 -out {}".format(gapFile.name, target, gapAlnFile.name)
    subprocess.call(remapCommand.split())
    alnFile = open(gapAlnFile.name)
    nMatch = 0
    matchLines = []
    matchStrands = []
    matchPositions = []
    matchValues = []
    matchPositions = []
    for alnLine in alnFile:
        alnVals = alnLine.split()
        chrom = alnVals[0]
        if (alnVals[3] == "0"):
            alnStart = int(alnVals[6])
            alnEnd  =  int(alnVals[7])
        else:
            alnStart = int(alnVals[8]) - int(alnVals[7])
            alnEnd = int(alnVals[8]) - int(alnVals[6])
        matchStrands.append(alnVals[3])
        seqLen = int(vals[4])

        if (alnEnd - alnStart  < 0.1 *seqLen):
            continue
        matchLines.append(alnLine)
        matchPositions.append((alnStart, alnEnd))
        matchValues.append(alnVals)
        nMatch += 1
        
    if (nMatch == 2 and matchStrands[0] != matchStrands[1]):
        sys.stdout.write("".join([ vals[3]+ " " + l for l in matchLines]))
        #
        # print the separation.
        print "{} {} {} {}".format(matchPositions[0][1], matchPositions[0][1],matchPositions[1][0], matchPositions[1][1])
        alnGap = min(abs(matchPositions[0][1] - matchPositions[1][0]), abs(matchPositions[1][0] - matchPositions[0][1]))
        
        print "aln gap: " + str(alnGap)
        command = "RepeatMasker " + gapFile.name
        
        devnull = open(os.devnull, 'w')
        subprocess.call(command.split(), stderr=devnull, stdout=devnull)
        outFile = open(gapFile.name + ".out")
        outFile.readline()
        outFile.readline()
        outFile.readline()
        repeatAnnotation=";".join( [ l.split()[9]+"/"+str(abs(int(l.split()[6]) - int(l.split()[5])))+"/"+l.split()[1] for l in outFile.readlines()])
        if (repeatAnnotation == ""):
            repeatAnnotation = "NONE"
        if (summaryFile is not None):
            chrom = chrom.split(":")[0]
            summaryFile.write(",".join(str(i) for i in [vals[3], chrom, min(matchPositions[0][1] - matchPositions[0][0], matchPositions[1][1] - matchPositions[1][0]), alnGap, repeatAnnotation]) + "\n")
        sys.stdout.write(repeatAnnotation + "\n")
        command=[ "/bin/rm" ] + [ a[0]+a[1] for a in zip([gapFile.name]*4, [".cat", ".masked", ".out", ".stderr"] ) ]
        subprocess.call(command)

        
    gapFile.close()
    if (args.deletetemp == True):
        subprocess.call(['/bin/rm', gapFile.name])
    else:
        print gapAlnFile.name

if (args.deletetemp == False):
    print refSeq.name
    print samFile.name
    print gapBedFile.name
    print gapFile.name
