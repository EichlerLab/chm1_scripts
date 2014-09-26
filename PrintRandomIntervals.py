#!/usr/bin/env python

import numpy as np
import sys
import Tools
import bisect

if (len(sys.argv) < 3):
    print "usage: PrintRandomIntervals fai nIntervals output"
    sys.exit(1)
    
fai = Tools.ReadFAIFile(sys.argv[1])

nIntv = int(sys.argv[2])
outFile = open(sys.argv[3], 'w')
cumsum = np.cumsum(np.asarray([c[0] for c in fai.values()]))

for i in range(nIntv):
    r = np.random.randint(0,cumsum[-1])
    ri = bisect.bisect_left(cumsum, r)
    chrom = fai.keys()[ri]
    pos   = np.random.randint(0,fai[chrom][0]-1)
    outFile.write("{}\t{}\t{}\n".format(chrom, pos, pos+1))
outFile.close()                  
    
