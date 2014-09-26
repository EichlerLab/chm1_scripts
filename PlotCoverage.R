#!/usr/bin/env Rscript

require(grDevices)
require(RColorBrewer)

BedOverlap <- function(bedTable, chrom, start, end) {
  chromi <- which(bedTable$V1 == chrom)
  newtab=bedTable[chromi,]
#  idx <- which((newtab$V2 <= start & newtab$V3 >= start) | (newtab$V2 <= end & newtab$V3 > end) | )
   idx <- which((newtab$V2 <= start & newtab$V3 >= start) | (newtab$V2 <= end & newtab$V3 > end )  | (newtab$V2 <= start & newtab$V3 >= end))
  return (newtab[idx,])
}

GetXBounds <- function(intv, refIntv) {
  return (c(max(refIntv[1], intv[1]), min(refIntv[2], intv[2])))
}


#
# Read specs
#
args <- commandArgs(trailingOnly=T)


if (length(args) != 8) {
  print("Usage: PlotCoverage.R coverageFile region segdup.bed assembly.bed chrom start end outputFile")
  quit("no")
}

region <- args[2]
chrom <- args[5]
start <- args[6]
end <- args[7]


#
# Read data
#
bed <- read.table(args[3],header=F)
asm <- read.table(args[4],header=F) 
cov <- read.table(args[1],header=F)


bed$V2=as.numeric(bed$V2)
bed$V3=as.numeric(bed$V3)
start=as.numeric(start)
end=as.numeric(end)


bedSubset <- BedOverlap(bed, chrom, start, end)
print(bedSubset)
nSubset = dim(bedSubset)[1]
p <- brewer.pal(min(11,max(3,nSubset)), "BrBG")

asmSubset <- BedOverlap(asm, chrom, start, end)
nAsm = dim(asmSubset)[1]
pasm <- brewer.pal(min(11, max(3,nAsm+1)), "RdYlGn")


pdf(args[8], width=4,height=5)


par(mar=c(3.2,3.1,1.5,2))
#layout(matrix(c(1,2,3), 3,1,byrow = TRUE), widths=c(3,3,3), heights=c(3,.75+.25*(nSubset ), .75+.25*(nAsm )))
layout(matrix(c(1,2,3), 3,1,byrow = TRUE), widths=c(3,3,3), heights=c(3,2,2))

plot(cov$V1,cov$V2,type='l',xlab="Position", ylab="Coverage",main=region , ylim=c(0,max(cov$V2)*1.10))

par(mar=c(4.5,3.1,.5,2))
plot(1,type="n", xlim=c(cov$V1[1], cov$V1[length(cov$V1)]), ylim=c(0,nSubset+1), axes=FALSE, xlab="Segmental duplication")
axis(1)
if (nSubset > 0) {
  for (i in seq(1,nSubset)) {
    xb <- GetXBounds(c(bedSubset$V2[i], bedSubset$V3[i]), c(cov$V1[1], cov$V1[length(cov$V1)]))
    rect(xb[1], .5+i-1, xb[2], 1.2+i-1, col=p[i%%length(p)])
    text(mean(xb),1+i-1.15, bedSubset$V4[i])
  }
}

par(mar=c(4.5,3.1,.5,2))
plot(1,type="n", xlim=c(cov$V1[1], cov$V1[length(cov$V1)]), ylim=c(0,nSubset+1), axes=FALSE, xlab="Assemblies")
axis(1)

if (nAsm > 0) {
  for (i in seq(1,nAsm)) {
    xb <- GetXBounds(c(asmSubset$V2[i], asmSubset$V3[i]), c(cov$V1[1], cov$V1[length(cov$V1)]))
    rect(xb[1], .5+i-1, xb[2], 1.2+i-1, col=pasm[i%%length(pasm)])
    text(mean(xb),1+i-1.15, asmSubset$V4[i])
  }
}




dev.off()


