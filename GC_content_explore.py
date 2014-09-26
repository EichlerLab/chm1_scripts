#!/usr/bin/env python
import pysam
import re
import matplotlib
import matplotlib.pyplot as plt
import numpy

from pbcore.io import CmpH5Reader
from pbcore.io import CmpH5Alignment

def IdentityFromCIGAR(cigar):
    nMatch = 0
    nIns   = 0
    nDel   = 0
    for cig in cigar:
        if (cig[0] == 0):
            nMatch += cig[1]
        elif(cig[0] == 1):
            nIns   += cig[1]
        elif(cig[0] == 2):
            nDel   += cig[1]
    denom = float(nMatch + nIns + nDel)
    return nMatch / denom

class  AlignmentSummary:
    def __init__(self, identity, length):
        self.identity = identity
        self.length   = length
        self.zmw      = 0
        self.sub      = 0

    
def SamToMap(samFileName, samMap):
    sf = pysam.Samfile( samFileName, "r" )
    for aln in sf.fetch():
        if (aln.rname == "*"):
            continue
        ident = IdentityFromCIGAR(aln.cigar)
        samMap[aln.qname] = AlignmentSummary(ident, aln.qlen)
        
def GetSubreadGC(subread):
    return (float(subread.basecalls().count('G') + subread.basecalls().count('C')) / len(subread.basecalls()))

def GetGC(read):
    maxLen = 0
    maxS   = 0
    for s in range(0,len(read.subreads)):
        l = len(read.subreads[s].basecalls())
        if (l > maxLen):
            maxLen = l
            maxS   = s
    return (float(read.subreads[maxS].basecalls().count('G') + read.subreads[maxS].basecalls().count("C"))) / len(read.subreads[maxS].basecalls())

        
        
#dh5 = "/net/eichler/vol20/projects/pacbio/backups/incoming/130625_MYD_eee_20kb_368/D01_1/Analysis_Results/m130626_034031_42134_c100534392550000001823079711101324_s1_p0.bas.h5"
#dsam = "/net/eichler/vol20/projects/pacbio/nobackups/results/130625_MYD_eee_20kb_368/D01_1/D.sam"

dh5 = "/mnt/pacbio/D01_1/Analysis_Results/m130626_034031_42134_c100534392550000001823079711101324_s1_p0.bas.h5"
dsam = "/mnt/pacbio_analysis/D01_1/D.sam"


from pbcore.io import BasH5Reader

dReader = BasH5Reader(dh5)

#
# key: 
#   rs  read score
#   rl  read length
#   mi  mapped identity 
#   ml  mapped length
#   m  indices of mapped reads
#   um  indices of unmapped reads
#   s  mapped subreads
#   us  unmapped subreads

class Count:
    def __init__(self):
        self.fields = ["rs", "rl", "mi", "ml", "m", "um", "s", "us"]
        self.data = { f: [] for f in self.fields }
        self.npdata = {}

    def ToNumpy(self):
        self.npdata = { f: numpy.array(self.data[f]) for f in self.fields }

def StoreMapped(fileNames, alnMap, stats):
    for fileName in fileNames:
        reader = BasH5Reader(fileName)
        for zmw in reader.sequencingZmws:
            for s in reader[zmw].subreads:
                stats.data["rs"].append(reader[zmw].readScore)
                stats.data["rl"].append(s.readEnd - s.readStart)
                if (s.readName in alnMap):
                    stats.data["m"].append(len(stats.data["rs"]) - 1 )
                    stats.data["ml"].append(alnMap[s.readName].length)
                    stats.data["mi"].append(alnMap[s.readName].identity)
                    stats.data["s"].append(s)
                else:
                    stats.data["um"].append(len(stats.data["rs"]) - 1 )
                    stats.data["ml"].append(0)
                    stats.data["mi"].append(0)
                    stats.data["us"].append(s)


dfn = ["/mnt/pacbio/D01_1/Analysis_Results/m130626_034031_42134_c100534392550000001823079711101324_s1_p0.bas.h5"]
dsam = "/mnt/pacbio_analysis/D01_1/D.sam"
dcmp = "/mnt/pacbio_analysis/D01_1/D.cmp.h5"

gfn = ["/mnt/pacbio/G01_1/Analysis_Results/m130626_103730_42134_c100534392550000001823079711101327_s1_p0.bas.h5","/mnt/pacbio/G01_1/Analysis_Results/m130626_103730_42134_c100534392550000001823079711101327_s2_p0.bas.h5"]
gsam = "/mnt/pacbio_analysis/G01_1/G.sam"

hfn = ["/mnt/pacbio/H01_1/Analysis_Results/m130626_125440_42134_c100534382550000001823079711101330_s1_p0.bas.h5","/mnt/pacbio/H01_1/Analysis_Results/m130626_125440_42134_c100534382550000001823079711101330_s2_p0.bas.h5"]
hsam = "/mnt/pacbio_analysis/H01_1/H.sam"

ffn = ["/mnt/pacbio/F01_1/Analysis_Results/m130626_081902_42134_c100534392550000001823079711101326_s1_p0.bas.h5","/mnt/pacbio/F01_1/Analysis_Results/m130626_081902_42134_c100534392550000001823079711101326_s2_p0.bas.h5"]
fsam = "/mnt/pacbio_analysis/F01_1/F.sam"

dStats = Count()
dh5Files = [dh5]
dSamMap = {}
SamToMap(dsam, dSamMap)
StoreMapped(dfn, dSamMap, dStats)
dStats.ToNumpy()

fStats = Count()
fSamMap = {}
SamToMap(fsam, fSamMap)
StoreMapped(ffn, fSamMap, fStats)
fStats.ToNumpy()

gStats = Count()
gSamMap = {}
SamToMap(gsam, gSamMap)
StoreMapped(gfn, gSamMap, gStats)
gStats.ToNumpy()

hStats = Count()
hSamMap = {}
SamToMap(hsam, hSamMap)
StoreMapped(hfn, hSamMap, hStats)
hStats.ToNumpy()

def ArrayHist(array, nbins=30):
    h = numpy.histogram(array, bins=nbins)
    return (h[1][0:-1], h[0])



def StatsHist(stats, dataset="rs", which="m", minValue=None):
    d = stats.npdata[dataset][stats.npdata[which]]
    if (minValue is not None):
        d = d[d > minValue]
    h = numpy.histogram(d, bins=30)
    return (h[1][0:-1], h[0])

dh = StatsHist(dStats, dataset="rs", which="m", minValue = 0.25)
fh = StatsHist(fStats, dataset="rs", which="m", minValue = 0.25)
duh =StatsHist(dStats, dataset="rs", which="um", minValue = 0.25)
fuh =StatsHist(fStats, dataset="rs", which="um", minValue = 0.25)

ax = plt.axes
plt.scatter(dh[0], dh[1], axes=ax)
plt.scatter(fh[0], fh[1], axes=ax, color="red")
plt.scatter(duh[0], duh[1], axes=ax, color="LightBlue")
plt.scatter(fuh[0], fuh[1], axes=ax, color="pink")
plt.show()

dCmpR = CmpH5Reader(dcmp)

mgc = numpy.array([GetSubreadGC(sr) for sr in gStats.npdata["s"]])
umgc = numpy.array([GetSubreadGC(sr) for sr in gStats.npdata["us"]])

dmgc = numpy.array([GetSubreadGC(sr) for sr in dStats.npdata["s"]])
dumgc = numpy.array([GetSubreadGC(sr) for sr in dStats.npdata["us"]])

hmgc = numpy.array([GetSubreadGC(sr) for sr in hStats.npdata["s"]])
humgc = numpy.array([GetSubreadGC(sr) for sr in hStats.npdata["us"]]
)

def GetLengths(subreads):
    return numpy.array([len(sr.basecalls()) for sr in subreads])
 
def IMean(array, indices):
    return np.mean(array[indices])

def LimitIndices(array, minValue = 0, maxValue=10000000):
    lowi = array > minValue
    highi = array < maxValue
    return lowi & highi



hl = GetLengths(hStats.npdata["s"])
hul = GetLengths(hStats.npdata["us"])

dl = GetLengths(dStats.npdata["s"])
dul = GetLengths(dStats.npdata["us"])

gmgch = ArrayHist(mgc)
gumgch = ArrayHist(umgc)

dmgch = ArrayHist(dmgc)
dumgch = ArrayHist(umgc)



ax1 = plt.subplot(121)
ax1.scatter(dl, dmgc, color="DarkRed", alpha=0.10)
ax1.scatter(hl, hmgc, color="DarkBlue", alpha=0.10)
ax2 = plt.subplot(122)
ax2.scatter(dul, dumgc, color="HotPink", alpha=0.10)
ax2.scatter(hul, humgc, color="DodgerBlue", alpha=0.10)
plt.show()


hi = LimitIndices(hmgc, 0.1, 0.7)
hui = LimitIndices(humgc, 0.1, 0.7)
di = LimitIndices(dmgc, 0.1, 0.7)
dui = LimitIndices(dumgc, 0.1, 0.7)

hi = hl > 1000
hui = hul > 1000
di =dl > 1000
dui = dul > 1000
ax1 = plt.subplot(121)
ax1.scatter(dl[di], dmgc[di], color="DarkRed", alpha=0.10)
ax1.scatter(hl[hi], hmgc[hi], color="DarkBlue", alpha=0.10)
ax2 = plt.subplot(122)
ax2.scatter(dul[dui], dumgc[dui], color="HotPink", alpha=0.10)
ax2.scatter(hul[hui], humgc[hui], color="DodgerBlue", alpha=0.10)
plt.show()

print numpy.mean(hmgc[hi])
print numpy.mean(humgc[hui])
print numpy.mean(dmgc[di])
print numpy.mean(dumgc[dui])

def GetGCContentByLength(lens, gc, nBins = 100):
    maxLength = np.max(lens)
    binSize   = maxLength/nBins
    gcBins    = [ [] for i in range(0,nBins)]
    
    for i in range(0,len(lens)):
        binIndex = min(int(lens[i]/binSize), nBins-1)
        gcBins[binIndex].append(gc[i])
    
    means = [ np.mean(gcBins[i]) if (len(gcBins[i]) > 0) else 0 for i in range(0,nBins) ]
    sds   = [ np.std(gcBins[i]) if (len(gcBins[i]) > 0) else 0 for i in range(0,nBins) ]
    x = [ binSize * i for i in range(0,nBins) ]
    return (x, np.array(means), np.array(sds))

(dx,dm,ds) = GetGCContentByLength(dl, dmgc)
(dux,dum,dus) = GetGCContentByLength(dul, dumgc)
(hx,hm,hs) = GetGCContentByLength(hl, hmgc)
(hux,hum,hus) = GetGCContentByLength(hul, humgc)

fig = plt.figure(figsize=(12,6))
ax1 = plt.subplot(121)
ax1.errorbar(hx,hm,yerr=hs, ecolor="DodgerBlue", color="blue")
ax1.errorbar(dx,dm,yerr=ds, ecolor="HotPink", color="red")
ax1.legend(("Dra1", "control"))
ax1.set(title="GC content of mapped reads by length")
ax1.axis([-1000,20000,0.1,0.7])
ax2 = plt.subplot(122)
ax2.errorbar(hux,hum,yerr=hus, ecolor="DodgerBlue", color="blue")
ax2.errorbar(dux,dum,yerr=dus, ecolor="HotPink", color="red")
ax2.set(title="GC content of unmapped reads by length")
ax2.axis([-1000,20000,0.1,0.7])
ax2.legend(("Dra1", "control"))
plt.show()

