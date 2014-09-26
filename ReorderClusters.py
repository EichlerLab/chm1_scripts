#!/usr/bin/env python

import sys
inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')

def Print(name, a, b):
    outFile.write(name)
    for val in a:
        outFile.write("\t" + val)
    for val in b:
        outFile.write("\t" + val)
    outFile.write("\n")
    
for line in inFile:
    vals = line.split()
    # the format is:
    # 0                                                                               1     2    3      4       5       6       7   8       9       10
    #m130216_080418_42134_c100465732550000001523050605101337_s1_p0/38995/5748_10266	-7740	0	chr1	16611	18883	-7803	1	chr1	16739	18961
    intv1 = (vals[1], vals[2], vals[3], vals[4], vals[5])
    intv2 = (vals[6], vals[7], vals[8], vals[9], vals[10].strip())
    if (int(vals[4]) < int(vals[9])):
        Print(vals[0], intv1, intv2)
    else:
        Print(vals[0], intv2, intv1)
