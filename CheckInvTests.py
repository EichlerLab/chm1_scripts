#!/usr/bin/env python
import sys
import argparse
import Tools

ap = argparse.ArgumentParser(description="Check for sequences that are better in the inverted orientation, using blasr m4 output.")
ap.add_argument("input", help="Alignments in m4 format.  This should be the alignment of sequences output  by DetectInversions.py")
ap.add_argument("--plot", help="Print commands to dotplot these reads. ", default=None)
ap.add_argument("--dir", help="Source dir of bas/bax.h5 files.", default=".")
args = ap.parse_args()

if (args.plot is not None):
    plotFile = open(args.plot, 'w')
    
inFile = open(args.input, 'r')
lines = inFile.readlines()

#m130928_232712_42213_c100518541910000001823079209281310_s1_p0/1331/0_13924/0_5772 chr8 -21290 85.7813 0 0 5772 5772 0 4905814 4911210 146364022 254 126251 -1835.77 -442 1

for i in range(len(lines)-1):
    v1 = lines[i].split()
    v2 = lines[i+1].split()
    t1 = v1[0].split('/')[0:3]
    t2 = v2[0].split('/')[0:3]
    s1 = int(v1[2])
    s2 = int(v2[2])
    if (t1 == t2):
        if (s2 < s1):
            diff = s2 - s1
            sys.stdout.write(str(diff) + " " + lines[i])
            sys.stdout.write(str(diff) + " " + lines[i+1])
            if (args.plot is not None):
                [barcode, zmw] = v1[0].split('/')[0:2]
                zmw = int(zmw)
                fileName = Tools.FindBasFile(barcode,zmw,args.dir)
                if (fileName is None):
                    print "Warning: could lot locatecontinue source data for " + barcode + " " + zmw
                    continue
                tStrand = int(v1[8])
                if (tStrand == 0):
                    start = int(v1[9])
                    end   = int(v1[10])
                else:
                    start = int(v1[11]) - int(v1[10])
                    end   = int(v1[11]) - int(v1[9])
                tName = v1[1]
                plotFile.write("~/projects/PacBioSequencing/scripts/DotPlot.py --query bas:{}:{} --target /var/tmp/mchaisso/ucsc.hg19.fasta --region {}:{}-{} --savefig DotPlots/{}_{}_{}.{}.{}.png --matches dot:11 \n".format(fileName, zmw, tName, start, end, tName, start, end, barcode, zmw))
                
                
