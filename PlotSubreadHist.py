#!/usr/bin/env python

import pysam
import sys


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import argparse


ap = argparse.ArgumentParser(description="Plot a histogram of subreads based on alignments.")
ap.add_argument("bam", help="Alignments of reads.")
ap.add_argument("--image", help="Write an image of the historgram here.", default=None)
ap.add_argument("--table", help="Write the table of the histogram here.", default=None)
ap.add_argument("--dname", help="Dataset name.", default=None)
opts = ap.parse_args()

samfile = pysam.Samfile( opts.bam, "rb" )
lengths = {}

for aln in samfile.fetch():
    idVals = aln.id.split('/')
    holeNumber = idVals[1]
    if (holeNumber not in lengths):
        lengths[holeNumber] = aln.qend - aln.qstart
    else:
        l = aln.qend - aln.qstart
        if (l < 
    lengths.append(aln.qend - aln.qstart)


fig = plt.figure()
ax  = plt.axes()
nplengths = np.asarray(lengths)
n, bins, patches = ax.hist(nplengths, 50, facecolor='green', alpha=0.75, log=True)
plt.figtext(.60,.75, "Median " + str(np.median(nplengths)))
plt.figtext(.60,.70, "95th " + str(np.percentile(nplengths, 95)))
plt.figtext(.60,.65, "Mean " + str(np.mean(nplengths)))
plt.figtext(.60,.60, "Total " + str(np.sum(nplengths)))
ax.set_xlabel("Subread lengths")
ax.set_ylabel("count")
if (opts.dname is not None):
    ax.set_title("Subread lengths " + opts.dname)
ax.grid(True)
plt.savefig(opts.image)



