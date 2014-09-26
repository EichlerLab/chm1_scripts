#!/usr/bin/env python

import sys
nameFile = open(sys.argv[1], 'r')

names = {}
for line in nameFile:
    vals = line.split()
    names[vals[0]] = True

first = True
for samFileName in sys.argv[2:]:
    samFile = open(samFileName, 'r')
    for line in samFile:
        if (first and len(line) > 0 and line[0] == '@'):
            sys.stdout.write(line)
        else:
            title = line.split()[0]
            if (title in names):
                sys.stdout.write(line)
    first = False
