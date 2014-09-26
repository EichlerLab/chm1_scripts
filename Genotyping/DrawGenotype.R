LoadFile <- function(filename) {
  return(read.table(filename, header=F))
}

#setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/Complex")

SaveMedianCount <- function(filename) {
  print(sprintf("reading %s\n", filename))
  t <- read.table(filename, header=F)
#  l <- unique(t$V1)
#  lidx <- seq(1,length(l))
#  lt <- length(t$V1)
#  i <- which(t$V1[1:(lt-1)] != t$V1[2:lt])
#
#  i <- c(0,i,lt)
#  print(i)  
#  idx <- sapply(seq(1,length(l)), function(i) which(t$V1 == l[i]))
#  return(sapply(seq(1,(length(i)-1)), function(j) mean(t$V3[i[j]:i[j+1]])))
  return(t$V2)
}

colorBins <- c(0,2,5, 10, 20, 40)

DrawGenotype <- function(genotype, col, row, pal, b=c(1,length(genotype)), intv=colorBins, verbose=0) {

  cols <- findInterval(genotype[b[1]:b[2]], intv)
  x0 <- seq(1,b[2]-b[1]+1) + col

  if (verbose != 0) {
    print("pal")
    print(pal[cols])
    print("genotype")
    print(genotype[b[1]:b[2]])
    print("intv")
    print(intv)
  }
  rect(x0, row, x0+1, row+1, col=pal[cols], border=NA)
  # redraw na

  naIdx <- which(is.na(genotype[b[1]:b[2]]))
  if (length(naIdx) > 0) {
    rect(col+naIdx, row, col+naIdx+1, row+1, col="red")
  }
}


DrawColorbar <- function(colors, x, y, len, hgt, labels, titleText="") {
  nColors <- length(colors)
  x0 <- rect(x+seq(0,nColors-1)*len, y, x+seq(1,nColors)*len,y+hgt, col=colors)
  par(srt=90)
  text(x+seq(0,nColors-1)*len+0.5, y+.3, labels, cex=0.25, pos=1)
  par(srt=0)
  text(x-1, y+hgt+2.5, labels=titleText, cex=0.5, pos=4)
}

require(getopt, quietly=T)
options <- matrix(c("input", "i", 2, "character",
                    "output", "o", 1, "character",
                    "queries", "q", 2, "character",
                    "roworder", "d", 2, "character",
                    "colorder", "e", 2, "character",
                    "chimp", "C", 2, "character",
                    "rows", "r", 2, "integer",
                    "columns", "c", 2, "integer"), byrow=T, ncol=4)

args <- getopt(options)

if (is.null(args$columns)) {
  columns <- 40
} else {
  columns <- args$columns
}

if (is.null(args$rows)) {
  rows <- 1
} else {
  rows <- args$rows
}

#
# Do some work to figure out what colors to paint the number of queries.
#



#
# Set up a coloring scheme for shading cells.
#

# number of matches per query
require(RColorBrewer, quietly=T)
fullPal <- brewer.pal(9,"Blues")
pal <- c("#FFFFFF", fullPal[5:9])

chimpPal <- c(0,brewer.pal(9,"Greens")[4:7])



#
# This has to run in the directory with all the count files in it.
#

# Read and compute the data to be plotted
fileTable <- read.table(args$input)
files <- paste(fileTable$V1)

names <- sapply(files, function(f) strsplit(f, '.', fixed=T)[[1]][1])
# plot the median count.

medians <- lapply(files, SaveMedianCount)

#
# store dimensions
#
nSamples <- length(medians)
nLoci <- length(medians[[1]])

if (is.null(args$roworder)) {
  rowOrder = seq(1,length(names))
} else {
  orderTable <- read.table(args$roworder)
  rowOrder <- orderTable$V1
  print("rowOrder")
  print(rowOrder)
}

if (is.null(args$colorder)) {
  colOrder = seq(1,nLoci)
} else {
  orderTable <- read.table(args$colorder)
  colOrder <- orderTable$V1
}



#
# Get the column names from the first genotype file.  They should all be the same.
#

f1 <- read.table(files[1])
queryCounts <- f1$V3
titles <- f1$V1
print("dim medians")
print(dim(medians))

medians <- medians[rowOrder]
names <- names[rowOrder]
print(sprintf("nloci %d", nLoci))
#medians <- medians[,colOrder]
#titles <- titles[colOrder]

if (is.null(args$queries) == F) {
  qt <- read.table(args$queries)
  queryOrder <- sapply(titles, function(i) match(i, qt$V1))
} else {
  queryOrder <- seq(1,length(titles))
}

# number of queries per cell

nQueryIntv <- 10^seq(0,ceiling(log10(max(queryCounts))))
nQueryPal <- c("#FFFFFF", brewer.pal(length(nQueryIntv)-1, "YlOrRd"))
countColors <- nQueryPal[findInterval(queryCounts, nQueryIntv)]
#
# If plotting if the sequence is present in chimp, do that here.
#
plotChimp <- F

genotypeOffset <- 0
if (is.null(args$chimp) == FALSE) {
#  args$chimp <- "ancestral_status/Complex.fasta.panTro4.bed"
#  chimp <- read.table("ancestral_status/Complex.fasta.panTro4.bed",header=F)
  chimp <- read.table(args$chimp, header=F)
  chimpIntv <- c(0,.25,.5,.75,1)
  chimpOrder <- sapply(titles, function(i) match(i,chimp$V1))
  plotChimp <- T
  genotypeOffset <- 1
  chimpOffset <- 1
}


#
# Determine the layout of the plots.
#


nStrips     <- ceiling(nLoci/columns)

nPages <- ceiling(nStrips / rows)
stripHeight <- 8
stripWidth  <- 6

labelTextColumns <- 7
labelPad <- 1.5
labelColumns <- labelTextColumns +  labelPad
layoutColumns <- columns + labelColumns
regionNumber <- 1

figHeight   <- stripHeight*rows 
print("fig etc")
print(c(rows))
pdf(args$output, width=stripWidth, height=figHeight)
outFileName <- sprintf("%s/%s",getwd(), args$output)
print(sprintf("Output in %s", outFileName))
par(mar=c(6,2,3,2))
par(xpd=TRUE)
colLabOffset <- -3
for (pageNo in seq(1,nPages)) {



  space <- 6
  plotStrips <- 1
  plotStrips <- min(nStrips - (rows * (pageNo-1)), rows)


#  dev.off()
#  x11(type="dbcairo", width=stripWidth, height=figHeight)


  maxY <- (nSamples+space)*plotStrips + 3
  plot(c(), xlim=c(0, layoutColumns), ylim=c(0,maxY), xlab="", ylab="", axes=F)
  
  for (i in seq(rows*(pageNo-1), rows*(pageNo-1)+plotStrips - 1 )) {
    b <- c(i*columns + 1, min((i+1)*columns, nLoci))

    yOffset <- (i - (rows*(pageNo-1)))*(nSamples+space) +2

#    text(2,at=c(seq(1,nSamples)+0.5 + yOffset),
    sapply(seq(1,nSamples), function(j) DrawGenotype(medians[[j]], labelColumns, j + yOffset + genotypeOffset, pal, b))
    text(1.5+colLabOffset,c(seq(1,nSamples)+yOffset+genotypeOffset+0.5),labels=(names), cex=0.25, adj=0, pos=4)
    
    sapply(seq(1,nSamples), function(j) DrawGenotype(queryCounts, labelColumns, yOffset, nQueryPal, b,intv=nQueryIntv))

    par(srt=90)

    text(labelColumns+1.25+seq(0,b[2]-b[1]), yOffset-6, sprintf("%s # %d", titles[b[1]:b[2]], queryOrder[b[1]:b[2]]), cex=0.5, pos=1)
    par(srt=0)
    if (plotChimp) {
      sapply(seq(1,nSamples), function(j) DrawGenotype(chimp$V4[chimpOrder], labelColumns, yOffset+chimpOffset, chimpPal, b,intv=chimpIntv))
      text(colLabOffset,yOffset+chimpOffset+.5, "Chimpanzee", cex=0.25, adj=0, pos=4)
    }

    text(colLabOffset,yOffset+0.5, labels="k-mer queries", cex=0.25, adj=0, pos=4)

    

  }
  colorBarOffset <- 1
  colorBarHeight <- 2
  DrawColorbar(pal, labelColumns+5, maxY+colorBarOffset, 1, colorBarHeight, paste(colorBins), "Read count")
  DrawColorbar(nQueryPal, labelColumns+5+8, maxY+colorBarOffset, 1, colorBarHeight, paste(nQueryIntv), "# k-mers")
  if (plotChimp) {
    DrawColorbar(chimpPal, labelColumns+5+16, maxY+colorBarOffset, 1, colorBarHeight, paste(chimpIntv), "panTro4 overlap")
  }
}
dev.off()
