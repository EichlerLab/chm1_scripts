#B!/usr/bin/env python
import sys
import Tools
import argparse
import subprocess
import time

ap=argparse.ArgumentParser(description="run a ton of dotplots")
ap.add_argument("mergedbed", help="Bed file with events.  The coordinates are of the event, and somewhere is a list of read names supporting the event, joined by ';'.")
ap.add_argument("genome", help="Genome file.")
ap.add_argument("outdir", help="Output dir")
ap.add_argument("readdir", help="Where to look for reads.")
ap.add_argument("--j", "-j", help="Number of threads", type=int, default=1)
ap.add_argument("--expand", help="Extend reference by this amount on either side", type=float, default=2.0)
ap.add_argument("--readcol", help="Get read names from this column (3)", type=int, default=3)
ap.add_argument("--titlecol", help="Read plot title from this column", default=None, type=int)
args = ap.parse_args()

inFile = open(args.mergedbed, 'r')
threads = []
for line in inFile:
    vals = line.split()
    reads = vals[args.readcol].split(";")
    readNames = []
    for read in reads:
        (barcode, zmw, coords) = Tools.ParseReadTitle(read)
        basFileName =  Tools.FindBasFile(barcode, int(zmw), args.readdir)
        if (basFileName is None):
            continue
        readNames.append("bas:{}:{}".format(basFileName, zmw))
    queries = " ".join(readNames)
    
    outName = args.outdir + "/" + vals[0]+"_"+vals[1]+"_"+vals[2] + ".pdf"

    if (args.titlecol is not None):
        title = "--title " + vals[args.titlecol]
    else:
        title = ""
    #
    # Check for an open thread
    #
    for t in range(len(threads)):
        if (threads[t].poll() != None):
            del threads[t]
            break
        if (len(threads) >= args.j):
            time.sleep(1)
    # reaching here means there are open thread positions
    gap = int(vals[2]) - int(vals[1])
    if (args.expand == 0):
        ext = int(gap*1.5) + 100
    else:
        ext = args.expand*gap
        
    command = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DotPlot.py --query {}  --target {} --matches dot:9 --savefig {} --region {}:{}-{} {}".format(queries, args.genome, outName, vals[0], int(int(vals[1]) - ext), int(int(vals[2])+ext), title)
    
    print "running " + command
    threads.append(subprocess.Popen(command.split()))
        


