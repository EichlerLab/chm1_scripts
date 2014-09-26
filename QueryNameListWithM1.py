#!/usr/bin/env python


import argparse
import sys
import bisect
import pickle
import Tools

ap = argparse.ArgumentParser(description="Query a name list using .m1 output")
ap.add_argument("m1", help="M1 input.")
ap.add_argument("table", help="Read Name table.")
ap.add_argument("--save", help="Save query here.", default=None)
ap.add_argument("--read", help="Read query dictionary.", default=False, action='store_true')
args = ap.parse_args()



if (args.read == True):
    dictFile = open(args.table, 'rb')
    table = pickle.load(dictFile)
else:
    tableFile = open(args.table, 'r')    
    table = dict( (line.split()[0], line.split()[1:]) for line in tableFile.readlines())
    
#m131004_073003_42213_c100572142530000001823103304021440_s1_p0/10/15165_18315	16	15217	18314	chr9	69878868	69881982
#table = [ v[0] for v in tableTxt ]
#table.sort()

if (args.save is not None):
    outFile = open(args.save, 'wb')
    pickle.dump(table, outFile)
    outFile.close()
    
alnFile = open(args.m1)


for line in alnFile:
    vals = line.split()
    query = vals[0]
    if (query not in table):
        print "Missing " + line.strip()
    else:
        #
        # look for truncated alignments
        aln = table[query]
        print str(aln)
        rStart = Tools.GetKV(aln, "XS:i:")
        rEnd   = Tools.GetKV(aln, "XE:i:")
        (bc,z,(hqs,hqe)) = Tools.ParseReadTitle(query)
        if (hqs is not None and hqe is not None and rStart is not None and rEnd is not None):
            if (hqe - rEnd  > 4000 or rStart - hqs > 4000):
                print "Truncated " + line.strip()
        else:
            print str((hqs, hqe, rStart, rEnd))
        
    
