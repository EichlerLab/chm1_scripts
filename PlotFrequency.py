#!/usr/bin/env python
#	Unrecognized backend string "asdf": valid strings are ['pgf', 'cairo', 'MacOSX', 'CocoaAgg', 'gdk', 'ps', 'GTKAgg', 'GTK', 'QtAgg', 'template', 'FltkAgg', 'emf', 'GTK3Cairo', 'GTK3Agg', 'WX', 'Qt4Agg', 'TkAgg', 'agg', 'svg', 'GTKCairo', 'WXAgg', 'pdf']
import matplotlib
#matplotlib.use('agg')
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

ap = argparse.ArgumentParser(description="Plot the nucleotide frequency of a pileup.")
ap.add_argument("--freq", help="Read the frequency from this file", required = True)
ap.add_argument("--out", help="Output file", required = True)

#args = ap.parse_args()


inFile = open(args.freq, 'rb')
freq = np.load("alignments.1.freq")

#freq = np.load(inFile)
margin = freq.sum(axis=0)


#plt.plot(freq[0])

margin +=1
frac = np.divide(freq, margin)
maxFrac = frac.max(axis=0)
