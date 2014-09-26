#!/usr/bin/env python

import sys

if (len(sys.argv) > 1):
    print "{}:{}-{}".format(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    while (sys.stdin):
        line = sys.stdin.readline()
        if (line == ""):
            sys.exit(0)
        else:
            line.strip()
            vals = line.split()
            print "{}:{}-{}".format(vals[0], vals[1], vals[2])
