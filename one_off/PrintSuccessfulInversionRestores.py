#!/usr/bin/env python
import sys

inFile = open(sys.argv[1])

prevLine = ""
prevTitle = ""
prevScore = 0
for line in inFile.readlines():
    title = '/'.join(line.split('/')[0:3])
    vals = line.split()
    score = int(vals[2])
    if (title == prevTitle and line.find("inverted") >= 0):
        
        if (score + 500 < prevScore):
            start = int(vals[9])
            end   = int(vals[10])
            if (vals[8] == "1"):
                start = int(vals[11]) - int(vals[10])
                end   = int(vals[11]) - int(vals[9])
            region = vals[1] +":" + str(start) +"-" +str(end) + "  " + vals[8]
            print vals[1] + "\t" + str(start) + "\t" + str(end) + "\t" + prevTitle + "\t" + title + "\t" + str(prevScore) + "\t" + str(score)
#            print prevLine.strip()
#            print line.strip()
#            print ""
    prevScore = score
    prevTitle = title
    prevLine  = line
