#!/usr/bin/env python


from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import argparse
import Tools
import sys
import subprocess
import tempfile

ap = argparse.ArgumentParser(description="Run dotplot of inversions.")

ap.add_argument("inv", help="Inversion file")
ap.add_argument("targets", help="Targets file")
ap.add_argument("genome", help="Genome file")
args = ap.parse_args()

#gHandle = open(args.genome, "r")
#genomeDict = SeqIO.to_dict(SeqIO.parse(gHandle, "fasta"))
#
gFile = file(args.genome, 'r')
gFai = Tools.ReadFAIFile(args.genome + '.fai')

tFile = file(args.targets, 'r')
tFai = Tools.ReadFAIFile(args.targets + '.fai')

#tHandle = open(args.targets, "r")
#targetDict = SeqIO.to_dict(SeqIO.parse(tHandle, "fasta"))


invFile = open(args.inv, 'r')
for line in invFile:
    vals = line.split()
    if (vals[4] != "interior"):
        continue
    chrom = vals[0]
    start = int(vals[1]) - 100
    end   = int(vals[2]) + 100

    invStart = int(vals[5])
    invEnd = int(vals[6])
    
    
    targetName = vals[3]
    tRegion = Tools.ParseRegionStr(targetName.split('|')[0])
    if (tRegion is None):
        print "ERROR with target: " + targetName
        continue
    tmpFileNames = []

    targetName = vals[3].split('/')[0]
    quiverName  = vals[3].split('|')[0] + "|quiver"
    splitValues = targetName.split('|')[1].split('_')
    splitStart = int(splitValues[1]) 
    splitEnd   = int(splitValues[2])
    quiverStart = splitStart + invStart - 100
    quiverEnd   = splitStart + invEnd   + 100
    genomeStart = tRegion[1]+splitStart
    genomeEnd   = tRegion[1]+splitEnd
    genomeStr = Tools.ExtractSeq((chrom, start, end), gFile, gFai)
    quiverStr = Tools.ExtractSeq((quiverName, quiverStart, quiverEnd), tFile, tFai)
    
    print "length of : " + targetName + " " + str(len(genomeStr))
    quiverName = targetName + " quiver"
    print "length of : " + quiverName + " " + str(len(genomeStr))
#    quiverStr = targetDict[targetName][tRegionStart:tRegionEnd]
    
    genomeFileName = tempfile.mktemp(suffix=".fasta", dir=".")
    quiverFileName = tempfile.mktemp(suffix=".fasta", dir=".")
    tmpFileNames.append(genomeFileName)
    tmpFileNames.append(quiverFileName)
    genomeFile = open(genomeFileName, 'w')
    quiverFile = open(quiverFileName, 'w')
    genomeName = "{}:{}-{}".format(chrom,start,end)
    quiverName = "{}:{}-{} quiver".format(chrom,start,end)
    SeqIO.write(SeqRecord.SeqRecord(seq=Seq.Seq(genomeStr), id=genomeName,name="",description=""), genomeFile, "fasta")
    SeqIO.write(SeqRecord.SeqRecord(seq=Seq.Seq(quiverStr), id=quiverName,name="",description=""), quiverFile, "fasta")
    genomeFile.close()
    quiverFile.close()
    command="/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DotPlot.py --query {} --target {} --savefig {}_{}_{}.pdf --matches dot:11".format(quiverFileName, genomeFileName, chrom,start,end)
    print "running " + command
    subprocess.call(command.split())
    for tmpfile in tmpFileNames:
        command = "/bin/rm -f " + tmpfile
        subprocess.call(command.split())
                
    
    
#chr10   102356490       102357412       chr10:100000000-110000000|quiver_2300000_2350000/0_50000        interior


#chr10:100000000-110000000|quiver_2300000_2350000/0_50000
