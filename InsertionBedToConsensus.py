#!/usr/bin/env python

import argparse
import sys
import os
import tempfile
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import subprocess
ap = argparse.ArgumentParser(description="Try and create consensi of insertion evetns based on a BED that has all inesrtion sequences support in it.")
ap.add_argument("bed", help="BED file. The 5th column (index 4) should be the support, and 8-8+col[4] are sequences.")
ap.add_argument("--path", help="Path to gcon.", default="/net/eichler/vol5/home/mchaisso/software/bin")
ap.add_argument("--minSupport", help="Only call consensus with at least this many reads.", default=3, type=int)
ap.add_argument("--out", help="Write out here.", default=None)

args = ap.parse_args()

if (args.out == None):
    outFile = sys.stdout
else:
    outFile = open(args.out, 'w')

bedFile = open(args.bed)
#fastaName = tempfile.mktemp(suffix=".fasta", dir=".")
#consensusName = tempfile.mktemp(suffix=".consensus.fasta", dir=".")
fastaName = "support.fasta"
consensusName = "consensus"
for line in bedFile:
    vals = line.split()
    nSupport = int(vals[4])
    index = 0
    firstSeq = vals[6].split(";")[0]
    lengths = [len(seq) for seq in vals[6].split(";")]
    lengths.sort()
    medianLength = lengths[len(lengths)/2]
    for seq in vals[6].split(";"):
        if (len(seq) == medianLength):
            medianSeq = seq
            break
    
    if (nSupport >= args.minSupport and len(firstSeq) > 200 ):
        seqFile = open(fastaName, 'w')

        for seq in vals[6].split(";"):
            title = "{}_{}_{}_{}".format(vals[0], vals[1], vals[1], index)
            index += 1
            SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq), id=title, name="",description=""),seqFile,  "fasta")
            
        seqFile.close()

        command = args.path + "/gcon.py d {} --min_cov 2 -o consensus ".format(fastaName, nSupport, consensusName)
        consFile = "consensus.fa"
        try:
            print "Running " + command
            subprocess.call(command.split())
        except:
            print "using default anyway"
            consensus = medianSeq
        if (os.path.exists(consFile)):
            inFile = open("consensus.fa")
            try:
                print "opening consensus!"
                seq = SeqIO.read(inFile, "fasta")
                consensus = seq.seq.tostring()


                if (consensus == ""):
                    consenus=medianSeq
            except:
                consensus=medianSeq
            command = "rm consensus.fa"
            subprocess.call(command.split())
        else:
            print consFile + " does not exist"
            consensus = medianSeq
    else:
        consensus = medianSeq
    outFile.write("\t".join(vals[0:6]) + "\t" + consensus + "\t" + "\t".join(vals[6:]) + "\n")
    consensus = "RESET!!!!"
    
cmd = "/bin/rm {}".format(fastaName)
subprocess.call(cmd.split())

cmd = "/bin/rm {}".format(consensusName)
subprocess.call(cmd.split())

if (args.out != sys.stdout):
    args.out.close()
