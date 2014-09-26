#!/usr/bin/env python

import pysam
import sys
import argparse
import Tools
import pdb
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Select the read extending most into a gap.")
ap.add_argument("bam", help="Indexed  file of alignments.")
ap.add_argument("region", help="Region to align", nargs="+")
ap.add_argument("--first", help="Write the longest extending read to this file.  If filedir is specified, this will be a .findex file suitable for input into writeHDFSubset.", default=None)
ap.add_argument("--last", help="If you're not first, you're last.  Print all reads overlapping the region boundary that are not the longest extension to this file.", default=None)
ap.add_argument("--side", help="Side to find overhang. Inside implies within a region, not overlapping", choices=['left', 'right', 'inside'], default='right')
ap.add_argument("--anchorOpposite", help="Make sure the non-overhang side is full length, by specifying a minimum number of bases that are missing.", type=int, default=None)
ap.add_argument("--filedir", help="Look for read files in this dirctory, and write output in .findex format: file<tab>zmw",default=None)
ap.add_argument("--seq", help="Print the read sequence from the bam and not the read name.")
ap.add_argument("--summary", help="Print a summary of the seed read overlap alignment to this file.", default=None)
ap.add_argument("--ref", help="Reference sequence. This is necessary when computing accuracy.", default=None)
ap.add_argument("--bestn", help="Print the 'n' best alignments (1)", type=int, default=1)
ap.add_argument("--mapq", help="Minimal mapping quality", type=int, default=30)
args = ap.parse_args()


    
region = Tools.FormatRegion(args.region)

    
if (region is None):
    print "malformatted region " + ' '.join(args.region)
    sys.exit(1)
    
samFile = pysam.Samfile(args.bam)

maxOverhang = 0
maxOverhangRead = None
overhangReads = {}
fai = None
refFile = None

if (args.summary is not None):
    summaryFile = open(args.summary, "w")

def GetClipping(cigar):
    frontClip = 0
    backClip = 0
    i = 0
    while (i < len(cigar) and (cigar[i][0] == 4 or cigar[i][0] == 5)):
        if (cigar[i][0] == 4):
            frontClip = cigar[i][1]
            break
        i +=1
    
    j = len(cigar) - 1
    while (j > i and (cigar[j][0] == 4 or cigar[j][0] == 5)):
        if (cigar[j][0] == 4):
            backClip = cigar[j][1]
        j-=1
    return (frontClip, backClip)

def GetTLen(cigar):
    tlen = 0
    for i in range(0,len(cigar)):
        if (cigar[i][0] == 0 or cigar[i][0] == 2):
            tlen += cigar[i][1]
    return tlen
    
def GetTagValue(tags, key):
    for t in tags:
        if (t[0] == key):
            return t[1]
    return None

refFile = None
if (args.ref is not None):
    refFile = open(args.ref)
    fai = Tools.ReadFAIFile(args.ref+".fai")


overhangs = []    
for aln in samFile.fetch(region[0], region[1], region[2]):
    print aln.qname
    if (aln.flag & 256 != 0):
        continue
    if (aln.mapq < args.mapq):
        continue
#    print aln.qname
    (barcode, zmw,coords) = Tools.ParseReadTitle(aln.qname)

    for t in aln.tags:
        if (t[0] == 'XS'):
#            print t
            readStart = t[1] - coords[0]
            xs = t[1]
        if (t[0] == 'XE'):
            readEnd   = t[1] - coords[0]
            xe = t[1]

    readLength = coords[1] - coords[0]

    readLen = xe - xs
    clip = GetClipping(aln.cigar)

    if (aln.flag & 16 != 0):
        strand = 1
    else:
        strand = 0
    
    readStart += clip[0]
    readEnd -= clip[1]
#    print clip
    anchorOpposite = 0

    if (args.side == 'right'):

 #       if (strand == 0):
#            print "strand {} overhang {} readlen {} readstart {} readend {} slack {} refStart {} refEnd {} ".format(strand, readLen - readEnd, readLen, readStart, readEnd, region[2]- (aln.pos + aln.tlen), aln.pos, aln.pos + aln.tlen)
        readOverhang = readLen - readEnd - (region[2]- (aln.pos + aln.tlen))
        anchorOpposite = clip[0]
#        else:
#            print "strand {} overhang {} readlen {} readstart {} readend {} slack {} refStart {} refEnd {}".format(strand, readStart, readLen, readStart, readEnd, region[2]- (aln.pos + aln.tlen), aln.pos, aln.pos + aln.tlen )
#            readOverhang = readStart - (region[2]- (aln.pos + aln.tlen))
#            anchorOpposite = clip[1]

    elif (args.side == 'left'):

#        if (strand == 0):
#            print "strand {} overhang {} readstart {} readend {} slack {}".format(strand, readStart, readEnd,aln.pos - region[1] )                    
        readOverhang = readStart - (aln.pos - region[1])
        anchorOpposite = clip[1]
#        else:
#            print "strand {} overhang {} readstart {} readend {} slack {}".format(strand, readLen-readEnd, readStart, readEnd,aln.pos - region[1] )                    
#            readOverhang = readLen - readEnd - (aln.pos - region[1])
#            anchorOpposite = clip[0]
    elif (args.side == 'inside'):
        if (strand == 0):
            readOverhang = readStart - (aln.pos - region[1]) + aln.pos + readLength - readStart - region[2]
        else:
            readOverhang = aln.pos + aln.rlen + readStart - region[2] + readEnd - (aln.pos - region[1])

    if (args.anchorOpposite is not None):
        if (anchorOpposite > args.anchorOpposite):
            # skip adding this read since it is truncated on the non-overhang side
            continue
            
    nt = 0
    c = aln.cigar
    for j in range(0,len(c)):
        if (c[j][0] == 0 or c[j][0] == 2):
            nt += c[j][1]

    if (refFile is not None):
        refseq = Tools.ExtractSeq((region[0], aln.pos, aln.pos+nt), refFile, fai)
        maxOverhangIdent = Tools.SAMToAccuracy(aln.cigar, aln.seq, refseq) * 100;
    else:
        maxOverhangIdent = 0
        
    overhangs.append([readOverhang, aln.qname, aln.seq, readStart, readEnd, region[0], aln.pos, aln.aend, aln.mapq, maxOverhangIdent, strand])


    if (readOverhang > 0):
        overhangReads[aln.qname] = readOverhang


overhangs.sort(reverse=True)

if (len(overhangs) == 0 or overhangs[0][0] < 0):
    print "No reads found in region."
    sys.exit(0)
    
def PrintExtensions(readDict, filedir, outFileName):
    if (outFileName != "/dev/stdout"):
        outFile = open(outFileName, 'w')
    else:
        outFile = sys.stdout
    printed = {}
    for readName, extend in readDict.iteritems():
        if (filedir is not None):
            res = Tools.ParseReadTitle(readName)
            base = res[0] + "/" + str(res[1])
            if (base in printed):
                continue
            printed[base] = True
            if (res is not None):
                (barcode, zmw, coord) = (res[0], res[1], res[2])
            else:
                print readName
                continue
            basFile = Tools.FindBasFile(barcode, zmw, filedir)
            if (basFile is not None):
                outFile.write(basFile + "\t" + str(zmw) + "\n")
            else:
                print "Could not find bas file for " + str((barcode, zmw, filedir))
        else:
            outFile.write(readName + "\t" + str(extend) + "\n")
    if (outFileName != "/dev/stdout"):
        outFile.close()

    
    # print the output to files.

    
if (args.seq is not None):
    seqFile = open(args.seq, 'w')
    print len(overhangs)
    for i in range(0, min(args.bestn, len(overhangs))):
        s = Seq.Seq(overhangs[i][2])
        r = SeqRecord.SeqRecord(s, id=overhangs[i][1], name="", description="")
        SeqIO.write(r, seqFile, "fasta")
    seqFile.close()


    
if (args.summary and args.seq is not None):
    for i in range(0,min(args.bestn, len(overhangs))):
        t = overhangs[i][0:2] + overhangs[i][3:]
        summaryFile.write("\t".join([str(j) for j in t]) + "\n")
        print t
    summaryFile.close()
                    

    
if (args.first is not None):
    reads = dict( [(maxOverhangRead, maxOverhang) ])
    PrintExtensions(reads, args.filedir, args.first)

else:
    print "------------------------------"
    if (len(overhangs) > 0):
        print str(overhangs[0][0]) + " " + overhangs[0][1]


if (args.last is not None):
    if (maxOverhangRead in overhangReads):
        del overhangReads[maxOverhangRead]
    print "after"
    PrintExtensions(overhangReads, args.filedir, args.last)

    
