#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
for line in inFile:
	vals = line.split()
	targetChrom = vals[3].split(":")[0]
	if (vals[0] == targetChrom):
		sys.stdout.write(line)

