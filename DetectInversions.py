#!/usr/bin/env python


import argparse
import sys
import Tools
import Bio.Seq as Seq
import Bio.SeqIO as SeqIO
import Bio.SeqRecord as SeqRecord
import pdb
ap = argparse.ArgumentParser(description="Find inversions in SAM files.")
ap.add_argument("--sam", help="Input sam files.", nargs="+", required=True)
ap.add_argument("--out", help="Output file, - = stdout", default="-")
ap.add_argument("--stride", help="Skip every N files.", default=1, type=int)
ap.add_argument("--start", help="Starting file.", default=0, type=int)
ap.add_argument("--tests", help="Write test seqeunces out to this file. For every test there is the forward aligned sequence, and the forward sequence with the inversion reversed.")
args = ap.parse_args()


            
if (args.out != "-"):
    outFile = open(args.out, 'w')
else:
    outFile = sys.stdout


def CreateTestSequences(aln0, aln1):
    # first make sequence forward direction
    (barcode, zmw, coords) = Tools.ParseReadTitle(aln0.title)
    seq0 = Seq.Seq(aln0.seq)
    qStart0 = aln0.qStart
    qEnd0   = aln0.qEnd
    
    if (aln0.flag & 16 != 0):
        seq0 = seq0.reverse_complement()

    invStart = aln1.qStart
    invEnd   = aln1.qEnd
    if (invStart < qStart0):
        return None
    invOffset = invStart - qStart0
    invLength = invEnd - invStart
    if (invOffset > len(seq0)):
        print "ERROR! " + str(invOffset) + " is too large " + str(len(seq0))
        
    revSeq = seq0[0:invOffset] + seq0[invOffset:(invOffset+invLength)].reverse_complement() + seq0[invOffset+invLength:]
    return (seq0, revSeq)
    
import pdb
def CheckForInversion(alns):
    # must have two alignments
    if (len(alns) < 2):
        return None
    # must have high mapqv
    if (alns[0].mapqv < 10):
        return None
    for j in range(1,len(alns)):
        # same chromosome
#        pdb.set_trace()
        if (alns[j].tName == alns[0].tName and
            ((alns[j].flag & 16) != (alns[0].flag & 16)) and
            alns[j].tStart > alns[0].tStart and
            alns[j].tEnd < alns[0].tEnd and
            alns[j].qStart > alns[0].qStart and
            alns[j].qEnd  < alns[0].qEnd):
            start = alns[j].tStart
            end   = alns[j].tEnd
            test = CreateTestSequences(alns[0], alns[j])
            return (alns[j].tName, start, end, "interior", test, alns[j].qStart, alns[j].qEnd)
    return None
    for j in range(1,len(alns)):
        if (alns[j].tName == alns[0].tName and # same chrom
            ((alns[j].flag & 16) != (alns[0].flag & 16)) and
            min(abs(alns[j].tStart - alns[0].tEnd), abs(alns[0].tStart - alns[j].tEnd)) < 5000):
            return (alns[j].tName, alns[j].tStart, alns[j].tEnd, "broken", None, alns[j].qStart, alns[j].qEnd)
    return None
    
if (args.tests is not None):
    tests = open(args.tests, 'w')

for i in range(args.start, len(args.sam), args.stride):
    samFileName = args.sam[i]
    sys.stderr.write(samFileName + "\n")
    sf = open(samFileName, 'r')
    # skip the header
    line = '@'
    prevTitle = ""
    while (len(line) > 0 and line[0] == '@'):
        line = sf.readline()
    # now parse alignments
    intervals = []
    while (len(line) > 0):
        vals = line.split()
        title = vals[0]
        #0      1     2      3     4      5      6    7        8
        aln = Tools.SAMEntry(line)
        if (aln.title == prevTitle or prevTitle == ""):
            intervals.append(aln)
        if (aln.title != prevTitle and prevTitle != ""):
#            if (len(intervals) > 1):
#                intervals[0].PrintIntervals(sys.stdout)
#                intervals[1].PrintIntervals(sys.stdout)
            coords = CheckForInversion(intervals)

            if (coords is not None):
                outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(coords[0], coords[1], coords[2], prevTitle, coords[3], coords[5], coords[6]))
            if (coords is not None and coords[4] is not None and args.tests is not None):
                rec = SeqRecord.SeqRecord(seq=coords[4][0], id=prevTitle, name="",description="")
                SeqIO.write(rec, tests, "fasta")
                rec2 = SeqRecord.SeqRecord(seq=coords[4][1], id=prevTitle + "/inverted", name="",description="")
                SeqIO.write(rec2, tests, "fasta")
#            print ""
                
            intervals = []
            intervals.append(aln)
        prevTitle = aln.title
        prevAln = aln
        line = sf.readline()
        
outFile.close()                
                
                
            
        
if (args.tests is not None):
    tests.close()
