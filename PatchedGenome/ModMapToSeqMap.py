#!/usr/bin/env python

import Tools
import sys
if (len(sys.argv) != 4):
    print "Usage ModMapToSeqMap.py modmap.bed genome.fasta seqmap.bed"
    print "  modmap is in the format:"
    print "  chr start end seq modchr modstart modend"
    print "  The seqmap is in the format: "
    print "  chr modstart modend MODseq"
    print "  The genome file should have an .fai index."
    sys.exit(1)
    
inFile = open(sys.argv[1], 'r')
genomeFileName = sys.argv[2]
outFile = open(sys.argv[3], 'w')

faiFile = Tools.ReadFAIFile(genomeFileName  + ".fai")
genomeFile = open(genomeFileName, 'r')

prevModStart = None

insPos = []
insSeq = []
import pdb


numBookEnds = 0

def PrintInsertion(chrom, modStart, modEnd, ins, seqs, genomeFile, faiFile, outFile):
    if (chrom == ""):
        return
    

    if (len(ins) == 0):
        return

    ins.insert(0, modStart)

    for i in range(1,len(ins)-1):
        lenPre = ins[i] - ins[i-1]
        prePos = ins[i-1]
        preSeq = Tools.ExtractSeq((chrom, prePos, prePos+lenPre), genomeFile, faiFile)
        global numBookEnds        
        numBookEnds+=1
        outFile.write(chrom + "\t" + str(prePos) + "\t" + str(prePos + lenPre)  + "\t" + preSeq + seqs[i-1] + "\t" + str(lenPre) + "\t" + str(len(seqs[i-1])) + "\t0" + "\n")

    i = len(ins) - 1
    lenPre = ins[i] - ins[i-1]
    prePos = ins[i-1]
    preSeq = Tools.ExtractSeq((chrom, prePos, prePos+lenPre), genomeFile, faiFile)

    sufSeq = Tools.ExtractSeq((chrom, prePos + lenPre, modEnd), genomeFile, faiFile)
    
    outFile.write(chrom + "\t" + str(prePos) + "\t" + str(prePos + lenPre)  + "\t" + preSeq + seqs[i-1] + sufSeq + "\t" + str(lenPre) + "\t" + str(len(seqs[i-1])) + "\t" + str(len(sufSeq)) + "\n")


insSeq = []
insPos = []
prevChrom = ""
prevStart = 0
prevEnd   = 0
doFinalPrint = True
for line in inFile:
    vals = line.split()
    modChrom = vals[4]
#    print str(vals)
    curStart = int(vals[5])
    curEnd   = int(vals[6])
    doPrint = False

    if (modChrom == "."):
        PrintInsertion(prevChrom, prevStart, prevEnd, insPos, insSeq, genomeFile, faiFile, outFile)
        PrintInsertion(vals[0], int(vals[1]), int(vals[1]), [int(vals[1])], [vals[3]], genomeFile, faiFile, outFile)
        insPos = []
        insSeq = []
        doFinalPrint = False
    else:
        if (curStart != prevStart):

            PrintInsertion(prevChrom, prevStart, prevEnd, insPos, insSeq, genomeFile, faiFile, outFile)
            insPos = []
            insSeq = []
        if (modChrom != vals[0]):
            print "not sure what to do with mismatching chromosomes, probably bad bed format."
            continue
        if (int(vals[5]) > int(vals[1])):
            print str(vals[0:3] + vals[5:])
            print "Not sure what to do in this case: genomic region after insertion start"
            continue
        prevStart = curStart
        prevEnd   = curEnd
        prevChrom = vals[0]
        insPos.append(int(vals[1]))
        insSeq.append(vals[3])

        doFinalPrint = True

if (doFinalPrint):        
    PrintInsertion(prevChrom, prevStart, prevEnd, insPos, insSeq, genomeFile, faiFile, outFile)        


sys.stderr.write("Printed {} book-ends\n".format(numBookEnds))


