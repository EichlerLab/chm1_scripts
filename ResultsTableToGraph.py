#!/usr/bin/env python

import sys

inFile = open(sys.argv[1], 'r')
# ditch the header
inFile.readline()

graph = {}

while(inFile):
    line1 = inFile.readline()
    line2 = inFile.readline()
    try:
        v1 = line1.split()[1]
        v2 = line2.split()[1]
    except IndexError:
        break
        
    #
    # Use lexicographic ordering of verties.
    #
    minV = min(v1, v2)
    maxV = max(v1, v2)
    if (minV not in graph):
        graph[minV] = {}
    if (maxV not in graph[minV]):
        graph[minV][maxV] = 0
    graph[minV][maxV] += 1


for src in graph.keys():
    for dest in graph[src].keys():
        print src + " " + str( graph[src][dest]) + " " + dest
        
            
    
