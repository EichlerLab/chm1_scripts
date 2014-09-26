#!/usr/bin/env python


import Tools

import sys

if (len(sys.argv) == 2):
    regionStr = sys.argv[1]
    regionBed = Tools.ParseRegionStr(regionStr)
    print "{}\t{}\t{}".format(regionBed[0], regionBed[1], regionBed[2])

else:
    while (sys.stdin):
        regionStr = sys.stdin.readline()
        if (regionStr == ""):
            sys.exit(0)
        regionStr.strip()
        regionBed = Tools.ParseRegionStr(regionStr)
        print "{}\t{}\t{}".format(regionBed[0], regionBed[1], regionBed[2])
        
        

