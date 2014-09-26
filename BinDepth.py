#!/usr/bin/env python

import sys
total = 0
n = 0
rgnLength = 0
prev = None
start = None
while True:
    l = sys.stdin.readline()
    if (l == ""):
        break
    v = l.split()
    v[1] = int(v[1])
    v[2] = int(v[2])
    if (prev is None):
        start = v[1]
    if (prev is not None and (prev[0] != v[0] or prev[1] + 1 != v[1] or (v[1]%100 == 0 and n > 0))):
        sys.stdout.write(prev[0] + "\t" + str(start) + "\t" + str(prev[1] + 1) + "\t{:2.2f}".format(total/n) + "\n")
        start = v[1]
        total = v[2]
        n = 1
    else:
        n += 1
        total += v[2]
    prev = v
                         
        
