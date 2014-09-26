LoadFile <- function(filename) {
  return(read.table(filename, header=F))
}



SaveMedianCount <- function(filename) {
  print(filename)
  t <- read.table(filename, header=F)
  l <- unique(t$V1)
  lidx <- seq(1,length(l))
  idx <- sapply(seq(1,length(l)), function(i) which(t$V1 == l[i]))
  return(sapply(idx, function(i) median(t$V3[i[1]])))
}

colorBins <- c(0,5, 10, 20, 40)

DrawGenotype <- function(genotype, col, row, pal, b=c(1,length(genotype)), intv=colorBins) {
  cols <- findInterval(genotype[b[1]:b[2]], intv)
  x0 <- seq(1,b[2]-b[1]) + col
  x1 <- seq(1,b[2]-b[1])+1 + col
  rect(x0, row, x1, row+1, col=pal[cols])
}


DrawColorbar <- function(colors, x, y, len, hgt, labels, titleText="") {
  nColors <- length(colors)
  x0 <- rect(x+seq(0,nColors-1)*len, y, x+seq(1,nColors)*len,y+hgt, col=colors)
  par(srt=90)
  text(x+seq(0,nColors-1)*len+0.5, y-.1, labels, cex=0.25, pos=1)
  par(srt=0)
  text(x, y+hgt+0.25, labels=titleText, cex=0.5, pos=4)
}

require(getopt, quietly=T)
options <- matrix(c("input", "i", 1, "character",
                    "query", "q", 1, "character",
                    "output", "o", 1, "character",
                    "chimp", "C", 2, "character",
                    "rows", "r", 2, "integer",
                    "columns", "c", 2, "integer"), byrow=T, ncol=4)



args <- getopt(options)

if (is.null(args$columns)) {
  columns <- 50
} else {
  columns <- args$columns
}

if (is.null(args$rows)) {
  rows <- 3
} else {
  rows <- args$rows
}

#
# Do some work to figure out what colors to paint the number of queries.
#

queries <- read.table(args$query, header=F)
queryCounts <- table(queries$V1)



#
# Set up a coloring scheme for shading cells.
#

# number of matches per query
require(RColorBrewer, quietly=T)
pal <- brewer.pal(5,"Blues")

chimpPal <- c("red", "green")

# number of queries per cell
nQueryPal <- brewer.pal(9, "YlOrRd")
nQueryIntv <- seq(0,max(queryCounts),max(queryCounts)/(9-1))
countColors <- nQueryPal[findInterval(queryCounts, nQueryIntv)]

#
# This has to run in the directory with all the count files in it.
#

# Read and compute the data to be plotted
files <- list.files(pattern="*.counts")

names <- sapply(files, function(f) strsplit(f, '.', fixed=T)[[1]][1])
# plot the median count.
print("saving medians")
medians     <- sapply(files, SaveMedianCount)


#
# If plotting if the sequence is present in chimp, do that here.
#
plotChimp <- F

genotypeOffset <- 0
if (is.null(args$chimp) == FALSE) {
  filename <- files[1]
  t <- read.table(filename, header=F)
  titles <- unique(t$V1)
  chimp <- read.table(args$chimp, header=F)
  chimpIntv <- c(0,1)
  chimpOrder <- sapply(titles, function(i) match(i,chimp$V1))
  plotChimp <- T
  genotypeOffset <- 1
  chimpOffset <- 1
}


#
# Determine the layout of the plots.
#
nStrips     <- ceiling(dim(medians)[1]/columns)

nPages <- ceiling(nStrips / rows)
stripHeight <- 3
stripWidth  <- 6

labelTextColumns <- 7
labelPad <- 1.5
labelColumns <- labelTextColumns +  labelPad
layoutColumns <- columns + labelColumns




for (pageNo in seq(1,nPages)) {
  print(sprintf("plotting %d", pageNo))
  outputFileName=paste(args$output, ".", sprintf("%d", pageNo), ".pdf", sep="")

  nMedians <- dim(medians)[1]

  nSamples <- dim(medians)[2]

  space <- 4

  plotStrips <- min(nStrips - (rows * (pageNo-1)), rows)
  print("plot strips")
  print(plotStrips)
  figHeight   <- stripHeight*plotStrips
  print(c(stripHeight, plotStrips, figHeight))
#  dev.off()
#  x11(type="dbcairo", width=stripWidth, height=figHeight)
  pdf(outputFileName,width=stripWidth, height=figHeight)
  
  maxY <- (nSamples+space)*plotStrips + 3
  plot(c(), xlim=c(0, layoutColumns), ylim=c(0,maxY), xlab="", ylab="", axes=F)
  
  nSamples <- dim(medians)[2]

  for (i in seq(rows*(pageNo-1), rows*(pageNo-1)+plotStrips - 1 )) {
    b <- c(i*columns + 1, min((i+1)*columns, dim(medians)[1]))
    par(mar=c(1,6,2,2))
    
    yOffset <- (i - (rows*(pageNo-1)))*(nSamples+space) +2

#    text(2,at=c(seq(1,nSamples)+0.5 + yOffset), 
    sapply(seq(1,nSamples), function(j) DrawGenotype(medians[,j], labelColumns, j + yOffset + genotypeOffset, pal, b))
    text(1.5,c(seq(1,nSamples)+yOffset+genotypeOffset+0.5),labels=(names), cex=0.25, adj=0, pos=4)
    sapply(seq(1,nSamples), function(j) DrawGenotype(queryCounts, labelColumns, yOffset, nQueryPal, b,intv=nQueryIntv))

    if (plotChimp) {
      sapply(seq(1,nSamples), function(j) DrawGenotype(chimp$V2[chimpOrder], labelColumns, yOffset+chimpOffset, chimpPal, b,intv=chimpIntv))
      text(1.5,yOffset+chimpOffset+.5, "Chimpanzee", cex=0.25, adj=0, pos=4)
    }

    text(1.5,yOffset+0.5, labels="k-mer queries", cex=0.25, adj=0, pos=4)

    

  }
  DrawColorbar(pal, labelColumns+5, maxY-2, 1, 1, paste(colorBins), "Read count")
  DrawColorbar(nQueryPal, labelColumns+5+10, maxY-2, 1,1, paste(nQueryIntv), "# k-mer queries")
  dev.off()
}
