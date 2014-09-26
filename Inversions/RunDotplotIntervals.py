#!/usr/bin/env python

import subprocess
import sys
import os
if (len(sys.argv) != 2):
    print "usage: RunDotplotIntervals.py input.bed"

inFile = open(sys.argv[1])

pbs = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts"
bam = "/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/chm1.pacbio.bam"
hg  = "/var/tmp/mchaisso/ucsc.hg19.fasta"
for line in inFile:
    vals = line.split()
    dirName = "_".join(vals)

    if (os.path.exists(dirName) is False):
        os.mkdir(dirName)

    region = "{}:{}-{}".format(vals[0], vals[1], vals[2])
    command = "{}/RegionToFasta.py {} {} --out {}/reads.fasta --filter /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/Inversions/MultipleBreaks/NearInvAlignments/reads.txt ".format(pbs, bam, region, dirName)
    print command
    # blocking call
    subprocess.call(command.split())
    delta = 10000
    targetRegion = "{}:{}-{}".format(vals[0], max(int(vals[1])-delta,0), int(vals[2])+delta)
    command = "samtools faidx {} {} > {}/target.fasta".format(hg, targetRegion, dirName)
    print command
    subprocess.call(command, shell=True)

    command = "{}/one_off/DotPlotQueriesToTarget.py {}/reads.fasta {}/target.fasta {}".format(pbs, dirName, dirName, dirName)
    subprocess.call(command.split())

    

    
    

