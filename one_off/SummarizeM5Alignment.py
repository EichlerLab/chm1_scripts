#!/usr/bin/env python
import math
import sys

inFile = open(sys.argv[1])

alnLine = inFile.readline()

vals = alnLine.split()
query = vals[16]
aln   = vals[17]
target= vals[18]

i = 0
hplen = 6
nPATIns = 0
nPATDel = 0

nDel = query.count('-')
nIns = target.count('-')
nMismatch = 0

while (i < len(query)):
    if (aln[i] == '*'):
        if (query[i] != '-' and target[i] != '-' and query[i] != target[i]):
            nMismatch += 1
            i += 1
            continue
        j = i + 1
        while (j < len(query) and aln[j] == '*'):
            j+=1
        ql = i
        qr = j
        tl = ql
        tr = qr
        while (ql > 0 and ((query[i-1] == 'T' or query[i-1] == 'A') and query[ql-1] == query[i-1] )):
            ql-= 1
        while (tl > 0 and ((target[i-1] == 'T' or target[i-1] == 'A') and target[tl-1] == target[i-1])):
            tl-= 1
        while (qr < len(query) and ((query[j] == 'T' or query[j] == 'A') and query[qr] == query[j])): 
            qr+= 1
        while (tr < len(target) and ((target[j] == 'T' or target[j] == 'A') and target[tr] == target[j])):
            tr+= 1
        if (query[i] == '-'):
            indel = 'del'
        else:
            indel = 'ins'
        if (i - ql > hplen or i - tl > hplen or qr - j > hplen  or tr - j > hplen):
            patlen = max(i - ql, i - tl , qr - j,  tr - j)
            motif = 'pAT'
            print indel + " " + motif + " " + str(j - i) + " " + str(patlen)
            if (indel == 'del'):
                nPATDel += j-i
            else:
                nPATIns += j-i
        i = j
    else:
        if (query[i] != target[i]):
            nMismatch += 1
        i += 1

print "summary: " + "npATdel: " + str(nPATDel) + " npATins: " + str(nPATIns)
print "mm: " + str(nMismatch)
print "length: " + str(len(target) - nDel)
print "total del: " + str(nDel) + " ins: "+  str(nIns)

print "phred: {:2.1f}".format(-10*math.log10((nMismatch + nDel + nIns) / float(len(target))))
