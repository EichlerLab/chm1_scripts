#!/usr/bin/env python
#chrM    4988    4989    m131012_092645_42215_c100597762550000001823105905221423_s1_p0/5607/0_2162       950     left

import sys

inFile = open(sys.argv[1])
prev = None
prevVals = ()

def Starts(a, b, maxDist):
    ast = int(a[1])
    aen = int(a[2])
    bst = int(b[1])
    ben = int(b[2])
    if (ast >= ben):
        d = ben - ast
        if (d < maxDist):
            return ((bst, ast), (ben, aen))
    elif (bst >= aen):
        if (bst - aen < maxDist):
            return ((ast, bst), (aen, ben))
    else:
        return None

prevVals = None
for line in inFile:
    curVals = line.split()
    if (prevVals is not None and prevVals[0] == curVals[0] and curVals[3] == prevVals[3]):
        starts = Starts(prevVals,curVals,100000)
        if (starts is not None):
            s =  str(starts[0][0])
            e =  str(starts[1][1])
            print prevVals[0] + "\t" + s + "\t" + e + "\t" + curVals[3] + "\t" + str(0) + "\t" + "+" + "\t" + s + "\t" + e + "\t" + "255,0,0" + "\t" + "2" + "\t" + ",".join([str(starts[0][1] - starts[0][0]), str(starts[1][1] - starts[1][0])]) + "\t" + "0," + str(starts[1][0] - starts[0][0]) 
    prevVals = curVals
            
            
