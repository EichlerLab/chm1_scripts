#!/usr/bin/env python

import subprocess
import sys
import datetime

cmph5File = sys.argv[1]
regionFile = open(sys.argv[2],'r')
outputDir  = sys.argv[3]
for region in regionFile:
    region = region.strip()
    commentPos = region.find("#")
    if (commentPos != -1):
        region = region[0:commentPos]
        region = region.strip()
    if (region == ''):
        continue
    outputFile = outputDir + "/" + region + ".fasta"
    outputFile = outputFile.replace(':',"_")
    outputFile = outputFile.replace('-',"_")
    command = "quiver -r /var/tmp/mchaisso/ucsc.hg19.fasta --referenceWindow {} -o {} {} ".format(region, outputFile, cmph5File)

    before = datetime.datetime.now()
    print command
    subprocess.call(command.split())
    after  = datetime.datetime.now()
    delta = after - before
    print str(delta)
    


