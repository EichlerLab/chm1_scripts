#!/usr/bin/env python

import sys
from Bio import Seq
from Bio.Alphabet import generic_dna, generic_protein

bedFile = open(sys.argv[1])
bedOut  = open(sys.argv[2], 'w')


def CountRepeat(repeat, seq):
    i = 0
    count = 0
    while (i < len(seq)):
        index = seq.find(repeat, i)

        if (index == -1):
            break
        else:
            i = index + len(repeat)
            count += 1
    return count

for line in bedFile:
    vals = line.split()
    seq  = vals[5]
    seq = seq.upper()
    repeat = vals[6]
    repeat = repeat.split(";")[0]

    vals[5] = seq
    i = 0
    l = len(seq)
    count = 0
    endRep = repeat.find(")")
    repeat = repeat[1:endRep]
    repeatRC = Seq.Seq(repeat, generic_dna).reverse_complement().tostring()
    forCount = CountRepeat(repeat, seq)
    revCount = CountRepeat(repeatRC, seq)
#    print repeat + " " + str(forCount) + " " + repeatRC + " "+  str(revCount)
    count = max(forCount, revCount)
    vals.append(str(count))
    vals.append("{:2.2f}".format( float(count*len(repeat))/len(seq)))
    bedOut.write("\t".join(vals) +"\n")
bedOut.close()
    

