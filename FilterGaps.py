#!/usr/bin/env python


import sys

nucMap           = [0]*256
nucMap[ord('A')] = 0
nucMap[ord('a')] = 0
nucMap[ord('C')] = 1
nucMap[ord('C')] = 1
nucMap[ord('G')] = 2
nucMap[ord('g')] = 2
nucMap[ord('T')] = 3
nucMap[ord('t')] = 3

inFile = open(sys.argv[1], 'r')
outFile =open(sys.argv[1], 'w')

for line in inFile:
    vals = line.split()
    
    #parse this
    #chr10   blasr   insertion       134714696,134714858     162     -       0       seq AAAAAAAAAAAAAACAAAAGAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAACCAAGAAAAAAAAAAAAAAACAGAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGTGGGGGGG  m130216_032125_42134_c100465732550000001523050605101335_s1_p0/7465/6126_10035   chr10:134714696-134714858
    seq = vals[8]
    counts = [0]*4
    for i in range(len(seq)):
        counts[nucMap[seq[i]]]+=1
    keep = True
    for i in range(len(counts)):
        if ( float(counts[i])/len(seq) > 0.9):
            keep = False
            break
    if (keep):
        outFile.write(line)
outFile.close()        
    
