#!/usr/bin/env python
import sys
import argparse
import Tools

ap = argparse.ArgumentParser(description="Print hard-stop coordinates, with a little wiggle-room, and the read names.")
ap.add_argument("sam", help="Input sam files.", nargs="+")
ap.add_argument("--out", help="Output file, 'stdout' defaults to stdout", default="stdout")
ap.add_argument("--window", help="Print this window around a hard-stop sequence.", default=50)
ap.add_argument("--delta", help="Minimum truncation to trigger a hard-stop", default=500)

args = ap.parse_args()
titlei = 0
flagi = 1
tnamei = 2
tposi  = 3
starti = 5
endi = 6
mapqvi = 4
tleni = 9

if (args.out == "stdout"):
    outFile = sys.stdout
else:
    outFile = open(args.out, "w")

index = 0    
for samFileName in args.sam:
    samFile = open(samFileName, 'r')
    sys.stderr.write("{}\t{}\t{}\n".format(index, len(args.sam), samFileName))
    index +=1
    for line in samFile:
        if (line[0] == "@"):
            continue
        else:
            entry = Tools.ParseSamLine(line)
            if (entry[mapqvi] < 30):
                continue
            order = "primary"
            if (entry[flagi] & 256 != 0):
                order = "secondary"
            (barcode, zmw, (rStart, rEnd)) = Tools.ParseReadTitle(entry[titlei])
            if (rEnd - entry[endi] > args.delta):
                pos = entry[tposi] + entry[tleni]
                outFile.write("{}\t{}\t{}\t{}\t{}\t{}\tR\t{}\t{}\n".format(entry[tnamei], max(pos-args.delta,0), pos+args.delta , entry[titlei], rEnd - entry[endi], order, entry[starti], entry[endi]))
            elif (entry[starti] - rStart > args.delta):
                pos = entry[tposi]
                outFile.write("{}\t{}\t{}\t{}\t{}\t{}\tL\t{}\t{}\n".format(entry[tnamei], max(pos-args.delta,0), pos+args.delta , entry[titlei], entry[starti], order, entry[starti], entry[endi]))
if (outFile != sys.stdout):
    outFile.close()
                            
                
