#!/usr/bin/env python


import argparse
import pysam
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt


ap = argparse.ArgumentParser(description="Summarize various aspects of alignments in SAM format.")
ap.add_argument("sam", help="Input.", nargs="+")
ap.add_argument("-F", help="Exclude these flages.", default=256)
ap.add_argument("-f", help="Include these flags.", default=0xFFFF)
ap.add_argument("-v", help="Verbose output.", action='store_true', default=False)
ap.add_argument("--plotInsert", help="Plot insert size distribution to this file", default=None)
ap.add_argument("--plotReadlen", help="Plot read length distribution to this file", default=None)
args = ap.parse_args()

nReads = 0
lengths = []
rLengths = []
subLengths = []
alnIndex = 0;

prl = {}
fl  = {}
for samFile in args.sam:
    if (samFile.find('.sam') >= 0):
        mode='r'
    else:
        # bam file
        mode='rb'
    f = pysam.Samfile(samFile,mode)
    sys.stderr.write(samFile + "\n")
    
    for aln in f.fetch():

        if (aln.flag & args.F != 0):
            continue
        
        xe = 0
        xs = 0
        zmwId = '/'.join(aln.qname.split('/')[0:2])
        tmp = aln.qname.split('/')[2].split('_')
        subreadStart = int(tmp[0])
        subreadEnd   = int(tmp[1])
        subLengths.append(subreadEnd - subreadStart)
        for t in aln.tags:
            if (t[0]=='XS'):
                xs = t[1]
            if (t[0] =='XE'):
                xe = t[1]
        if (xe != 0):
            rLengths.append(xe - xs)
        if (args.v == True):
            print "xlen: " +str(xe - xs)
            print "len: " + str(aln.tlen)
            print "flag: " + str(aln.flag)
            print aln
        lengths.append(aln.tlen)
        if (zmwId not in prl):
            prl[zmwId] = 0
        prl[zmwId] += aln.tlen
        if (zmwId not in fl):
            fl[zmwId] = 0
        if (fl[zmwId] < aln.tlen):
            fl[zmwId] = aln.tlen
            
        if (alnIndex % 10000 == 0 and alnIndex > 0):
            sys.stderr.write("processed {}\n".format(alnIndex))
        alnIndex += 1
    f.close()


polyLengths= np.asarray(prl.values())
insertLengths = np.asarray(fl.values())
polyLengths.sort()
insertLengths.sort()    
npLengths = np.asarray(lengths)
nprLengths = np.asarray(rLengths)
npLengths.sort()
nprLengths.sort()
npsLengths = np.asarray(subLengths)
npsLengths.sort()

print "#Files {}".format(len(args.sam))
print "#Alignmetns {}".format(len(npLengths))
print "#Read Alignments {}".format(len(nprLengths))
print "#bases {}".format(sum(npLengths))
titles = ["RefAlignLength", "ReadAlignLength", "SubreadLength", "Proccessed", "Insert"]
counts = [len(npLengths), len(nprLengths), len(npsLengths), len(polyLengths), len(insertLengths)]
means = [npLengths.mean(), nprLengths.mean(), npsLengths.mean(), polyLengths.mean(), insertLengths.mean()]
medians = [npLengths[len(npLengths)/2], nprLengths[len(nprLengths)/2], npsLengths[len(npsLengths)/2], polyLengths[len(polyLengths)/2], insertLengths[len(insertLengths)/2]]
tops = [npLengths[int(len(npLengths)*0.95)], nprLengths[int(len(nprLengths)*0.95)], npsLengths[int(len(npsLengths)*0.95)], polyLengths[int(len(polyLengths)*.95)], insertLengths[int(len(insertLengths)*0.95)]]
total = [npLengths.sum(), nprLengths.sum(), npsLengths.sum(), polyLengths.sum(), insertLengths.sum()]
print "Type \t" + "\t".join(titles) + "\n"
print "N:    \t" + "\t".join([str(m) for m in counts])
print "Mean:  \t " + "\t".join(["{:.2f}".format(m) for m in means])
print "Median: \t" + "\t".join([str(m) for m in medians])
print "95th:   \t" + "\t".join([str(i) for i in tops])
print "Total:  \t" + "\t".join([str(s) for s in total])


