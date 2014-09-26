#!/usr/bin/env python

from Bio import SeqIO 
from Bio import SeqRecord
from Bio import Seq
import subprocess
import Tools
import sys
import time
import tempfile
import os
if (len(sys.argv) != 4):
    print "USage RunDotPlot.py reads.multifasta.fasta target.fasta outdir"
    sys.exit(1)
    
readsName  = sys.argv[1]
targetName = sys.argv[2]
outdir  = sys.argv[3]


reads = open(readsName)
j=10
threads = []
targetFile = open(targetName)
target = SeqIO.read(targetFile, "fasta")
tmpQueryNames = []

for query in SeqIO.parse(reads, "fasta") :
    tmpQueryName = tempfile.mktemp(suffix=".query.fasta", dir=".")
    tmpQueryNames.append(tmpQueryName)
    oneRead = open(tmpQueryName, 'w')
    SeqIO.write(query, oneRead, "fasta")

    oneRead.close()
    queryVals = query.id.split("_")
    readName = query.id.replace("/","_")
    dotPlotName = outdir + "/" + readName + ".png"
    dotPlotName = dotPlotName.replace(":","_")
    (chrom, start, end) = Tools.ParseRegionStr(target.id)
    xStart = start
    command="/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DotPlot.py --query {} --target {} --savefig {}  --matches dot:13 --thin  --xstart {} --ylabel {} --nolegend ".format(tmpQueryName, targetName, dotPlotName, xStart, query.id, queryVals[0])
    print command

    while (len(threads) >= j):
        for t in range(0,len(threads)):
            if (threads[t].poll() != None):
                del threads[t]
                del tmpQueryNames[t]
            break
        if (len(threads) >= j):
            time.sleep(1)

    threads.append(subprocess.Popen(command.split()))


