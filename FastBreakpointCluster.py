#!/usr/bin/env python


import sys

inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')


allClusters = {}
delta = 200

class Interval:
    def __init__(self, chr1, str1, sta1, chr2, end1, str2, sta2, end2):
        self.intv = [(chr1, str1, sta1, end1), (chr2, str2, sta2, end2)]
    
    def Overlaps(self, str1, sta1, end1, str2, sta2, end2):
        if (self.intv[0][0] != chr1 or
            self.intv[1][0] != chr2):
            return False
        
        if (self.intv[0][1] == str1 and
            self.intv[1][1] == str2 and
            ((self.intv[0][2] <= sta1 and self.intv[0][3] >= end1) or (abs(self.intv[0][2] - sta1) < delta) or (abs(self.intv[0][3] - end1) < delta)) and
            ((self.intv[1][2] <= sta2 and self.intv[1][3] >= end2) or (abs(self.intv[1][2] - sta2) < delta) or (abs(self.intv[1][3] - end2) < delta))):
            return true
    def Overlaps(self, intv):
        return self.Overlaps(intv[0][0], intv[0][1], intv[0][2], intv[0][3], intv[1][0], intv[1][1], intv[1][2], intv[1][3])

    
    def Merge(self, str1, sta1, end1, str2, sta2, end2):
        if (str1 < self.intv[0][2]):
            self.intv[0][2] = str1;
        if (end1 > self.intv[0][3]):
            self.intv[0][3] = end1
        if (sta2 < self.intv[1][2]):
            self.intv[1][2] = sta2
        if (end2 > self.intv[1][3]):
            self.intv[1][3] = end2

    def Merge(self, intv):
        self.Merge(intv[0][1], intv[0][2], intv[0][3], intv[1][1], intv[1][2], intv[1][3])
        
lines = inFile.readlines()
    # 0                                                                               1     2    3      4       5       6       7   8       9       10
    #m130216_080418_42134_c100465732550000001523050605101337_s1_p0/38995/5748_10266	-7740	0	chr1	16611	18883	-7803	1	chr1	16739	18961
vals = [line.split() for line in lines]
intervals = [Interval(v[3], int(v[2]), int(v[4]), int(v[5]), v[8], int(v[7]), int(v[9]), int(v[10])) for val in vals]
for i in range(0,len(lines)):
    j = i + 1;
    

    vals = line.split()
    lscore = int(vals[1])
    lstrand = int(vals[2])
    lchr = vals[3]
    lstart = int(vals[4])
    lend   = int(vals[5])
    rscore = int(vals[6])
    rstrand = int(vals[7])
    rchr = vals[8]
    rstart = int(vals[9])
    rend   = int(vals[10])
    if (lchr not in allClusters):
        allClusters[lchr] = []

    cluster = allClusters[lchr]
    
    
