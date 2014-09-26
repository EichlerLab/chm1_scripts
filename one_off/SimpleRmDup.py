#!/usr/bin/env python

import sys
prev = [0,0,0]
for line in sys.stdin:
    vals = line.split()
    if (vals[0] == prev[0] and vals[1] == prev[1] and vals[2] == prev[2]):
        continue
    sys.stdout.write(line)
    prev = vals
