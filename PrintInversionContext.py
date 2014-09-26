#!/usr/bin/env python

import argparse
import sys
import Tools
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="Take an inversion bed file, and print the sequences in inverted format with some flanking part of the genome.")
ap.add_argument("bed", help="Input BED file.  This should be in the insertion bed format.")
ap.add_argument("genome", help="Input genome with a .fai file.")
ap.add_argument("contexts", help="Output file.")
ap.add_argument("--window", help="Amount to store on sides.", type=int, default=1000)
ap.add_argument("--bponly", help="Print only the breakpoints, value specifies amount of sequence about each breakpoint to print.", type=int, default=None)
ap.add_argument("--noreverse", help="Assume the breakpoints are from assembled contigs and there is no need to reverse the alignment.", default=False, action='store_true')

args = ap.parse_args()

bedFile = open(args.bed)
contextFile = open(args.contexts, 'w')

fai = Tools.ReadFAIFile(args.genome+".fai")

genomeFile = open(args.genome)


for line in bedFile:
    vals  = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end = int(vals[2])
    
    prefixCoords = (max(0, start - args.window), start)
    suffixCoords = (end, min(fai[chrom][0], end+args.window))

    invSeq = Seq.Seq(Tools.ExtractSeq((chrom, start, end), genomeFile, fai))
    prefix = Tools.ExtractSeq((chrom, prefixCoords[0],prefixCoords[1]), genomeFile, fai)
    suffix = Tools.ExtractSeq((chrom, suffixCoords[0],suffixCoords[1]), genomeFile, fai)
    if (args.noreverse == False):
        invSeqStr =  invSeq.reverse_complement().tostring()
    else:
        invSeqStr =  invSeq.tostring()
        
    context = prefix + invSeqStr + suffix
    seqname = "/".join(vals[0:3])
    if (args.bponly is None):
	SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(context), id=seqname, name="",description=""), contextFile, "fasta")
    else:
	bp1 = len(prefix)
	bp2 = len(prefix) + len(invSeqStr)
	bp1Seq = context[bp1-args.bponly:bp1+args.bponly]
	bp2Seq = context[bp2-args.bponly:bp2+args.bponly]
	bp1Name = seqname + "_1"
	bp2Name = seqname + "_2"
	SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(bp1Seq), id=bp1Name, name="",description=""), contextFile, "fasta")
	SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(bp2Seq), id=bp2Name, name="",description=""), contextFile, "fasta")
    
contextFile.close()
    
