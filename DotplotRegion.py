#!/usr/bin/env python

import sys
import subprocess
import Tools
import os

if (len(sys.argv) > 2):
    region = Tools.FormatRegion(sys.argv[1:])
else:
    region = Tools.FormatRegion(sys.argv[1])
s = region[1] - 2000
e = region[2] + 2000

dirname = "{}_{}".format(region[0], str(region[1])) 
os.mkdir("plots/"+dirname)

command =  "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DotPlot.py --query region:{}:{}-{}  /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/chm1.pacbio.bam  --region {}:{}-{} --target /var/tmp/mchaisso/ucsc.hg19.fasta --savefig plots/{}/{}.png --maxq 12 --matches dot:11   --tappend /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/MEI/LINE_L1.fasta ".format(region[0], region[1], region[2], region[0], s, e, dirname, dirname)
subprocess.call(command.split())


