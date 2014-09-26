#!/usr/bin/env python

import argparse
import subprocess
import tempfile
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

ap = argparse.ArgumentParser(description="Given a gapped bed file, repeat mask the sequences in it, and create a BED file with the sequence as repeat masked.")


ap.add_argument("bedin", help="Input bed file")
ap.add_argument("bedout", help="Output bed file.")
ap.add_argument("--tmpdir", help="Working directory", default=None)


args = ap.parse_args()


if (args.tmpdir is None):
    args.tmpdir = tempfile.mktemp()
    
cmd = "mkdir -p {}".format(args.tmpdir)
subprocess.call(cmd.split())


#
# Print the fasta sequences
#
tmpSeqFileName = tempfile.mktemp(suffix=".tomask.fasta", dir=args.tmpdir)
source = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts"

command = "{}/GapBedToFasta.py {} {}".format(source, args.bedin, tmpSeqFileName)
print command
subprocess.call(command.split())

command = "RepeatMasker -xsmall -dir {} -pa 8 {}".format(args.tmpdir, tmpSeqFileName)
print command
subprocess.call(command.split())

# now read the masked file
maskedSequences = open(tmpSeqFileName+".masked")
maskedDict = SeqIO.to_dict(SeqIO.parse(maskedSequences, "fasta"))

inBedFile = open(args.bedin, 'r')
outBedFile = open(args.bedout, 'w')

for line in inBedFile:
    vals = line.split()
    title = '/'.join(vals[0:3])
    if (title in maskedDict):
        vals[5] = maskedDict[title].seq.tostring()
    outBedFile.write('\t'.join(vals) +"\n")

inBedFile.close()
outBedFile.close()

command = "/bin/rm -rf {} ".format(args.tmpdir)
subprocess.call(command.split())


    

