#!/usr/bin/env python

import sys

repeat = sys.argv[1]
seq    = sys.argv[2]

repeat = repeat.upper()
seq = seq.upper()

count = 0
i = 0
l = len(seq)
while (i < len(seq)):
    index = seq.find(repeat, i)
    if (index == -1):
        break
    else:
        i = index + len(repeat)
        count += 1


print str(count) +  "\t" + str(count*len(repeat)) + "\t" + str(len(seq))
sys.exit(0)

