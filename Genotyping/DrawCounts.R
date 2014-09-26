require(getopt, quietly=T)
options <- matrix(c("input", "i", 2, "character",
                    "output", "o", 2, "character",
                    "norm", "n", 1, "character",
                    "nrows", "r", 1, "integer",
                    "roi", "e", 1, "character",
                    "labels", "l", 1, "character"
                    ), byrow=T, ncol=4)

args <- getopt(options)

require(RColorBrewer, quietly=T)



mav <- function(x,n=50){filter(x,rep(1/n,n), sides=2)}

SaveCount <- function(filename, nr=-1) {
  print(sprintf("reading %s\n", filename))
  t <- read.table(filename, header=F,nrows=nr, colClasses=c("character", "NULL", "integer", "integer"))
  n <- names(table(t$V1))
  return(sapply(n, function(i) t$V3[which(t$V1 == i)]))
}

SavePositions <- function(filename, nr=-1) {
  print(sprintf("reading %s\n", filename))
  t <- read.table(filename, header=F,nrows=nr, colClasses=c("character", "NULL", "integer", "integer"))
  n <- names(table(t$V1))
  return(sapply(n, function(i) t$V4[which(t$V1 == i)]))
}





# testing code
#setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/STR/by_position")
#setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/Complex")


if (is.null(args$nrows)) {
  nRows = -1
} else {
  nRows = args$nrows
}
#nRows <- 10000
print(sprintf("Reading input table: '%s'", args$input))
files <- read.table(args$input)
print(files)
#files <- read.table("counts.test.txt")
#nRows <- 2178409
files <- paste(files$V1)

roiTable <- c()
#roiTable <- read.table("STR.fasta.period")
if (is.null(args$roi) == F) {
  roiTable <- read.table(args$roi)
  head(roiTable)
}

if (is.null(args$labels)) {
  labels <- c("WEA", "EA", "SA", "AFR", "ADM", "CHM1")
} else {
  
  lt <- read.table(args$labels)
  labels <- paste(lt$V1)
  print("labels")
  print(labels)

}

getBase <- function(i) {
  b=unlist(strsplit(i,"[.]"))
  return(b[1])
}


positions <- SavePositions(files[1], nr=nRows)
allCounts <- lapply(files, SaveCount, nr=nRows)

n <- names(allCounts[[1]])

if (is.null(args$norm)) {
  biasTable <- read.table("/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/CHM1Sequencing/Genotyping/control_regions/RegionNormalization.norm.txt")
} else {
  biasTable <- read.table(args$norm)
}


fileBase <- sapply(files, getBase)
#fileBase <- "CHM1.UW.counts"
bias <- sapply(fileBase, function(f) biasTable$V2[grep(f, biasTable$V1)])
bias = as.numeric(bias)
# flatline anything that is not matched.
bias[which(is.na(bias))] = 1


IndexedMean <- function(counts, i, name, indices) {
  aboveZero <- which(counts[[i]][[name]] > 0)
  l = length(aboveZero)
  if (l == 0) {
    return(0)
  } else {
    return(mean(counts[[i]][[name]][aboveZero]))
  }
}


GroupMean <- function(counts, name, indices, index) {
  print("indicex")
  print(indices)
  print("index")
  print(index)
  
  pi <- which(unlist(indices) == index)
  
  return(sapply(pi, function(ii) IndexedMean(counts, ii, name, pi)))
}

PlotBoxes <- function(counts, name, indices, maxIndex=5, colors=c("black"), legLabels=c(""), legPal=c("black")) {
  vals <- lapply(seq(1,maxIndex), function(i) GroupMean(counts, name, indices, i))
  print(vals)
  boxplot(vals, xlab=name, col=legPal)
  legend("topright", legend=legLabels, col=legPal, pch=15)
}


DrawCountPoints <- function(counts, positions, color, window) {
  windowSize <- min(window, length(counts))
  movingAverage <- mav(counts, windowSize)
  indices <- which(is.na(movingAverage) == F)
  points(positions[indices], movingAverage[indices], col=color, type='l')
#  points(positions, counts, col=color, type='l')
#  points(positions[indices], movingAverage[indices], col=color)
}

PlotCount <- function(counts, positions, name, window=10, colors=c("black"), legLabels=c(""), legPal=c("black"), rgn=c(0,0)) {
  nSeq <- length(counts)

  print(sprintf("plotting %s",name))
  nCounts <- length(counts)
  gm <- max(sapply(seq(1,nCounts), function(i) max(na.omit(counts[[i]][[name]]))))
  print(gm)
  print(sprintf("nelem %d", nCounts))
  plot(c(), ylim=c(0,gm), xlim=c(0, max(positions)), xlab="Position in insert", ylab="Normalized read count",main=name)

  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else{  DrawCountPoints(counts[[i]][[name]], positions, colors[i], window); }})
  print(rgn)
  if (length(rgn) == 2 && !is.na(rgn[1]) && !is.na( rgn[2]) ) {
    segments(rgn[1], 0, rgn[2], 0, lwd=2) }
  legend("topleft", legend=legLabels, col=legPal, lty=1)

}

PlotLog2 <- function(counts, positions, name, window=100, colors=c("black"), legLabels=c(""), legPal=c("black")) {
  print(sprintf("log plotting %s\n",name))
  nCounts <- length(counts)
  means <- sapply(seq(1,nCounts), function(i) mean(counts[[i]][[name]]))
  medianI <- order(means)[floor(length(means)/2)]
#  i <- 3
#  counts[[i]][[name]]/counts[[medianI]][[name]]
  
  logVals <-  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else { vals <-log(counts[[i]][[name]]/counts[[medianI]][[name]],2);  vals[which(is.na(vals))] = 0;  vals[which(is.infinite(vals))] = 0; return(vals)  }})
  
  plot(c(), ylim=c(min(sapply(logVals, min)), max(na.omit(sapply(logVals, max)))), xlim=c(0, max(positions)), xlab="Position in insert", ylab="Read count",main=name)
  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else { vals <-log(counts[[i]][[name]]/counts[[medianI]][[name]],2);  vals[which(is.na(vals))] = 1; DrawCountPoints(vals , positions, colors[i], window); }})

  legend("topleft", legend=legLabels, col=legPal, lty=1)

}



nCounts <- length(allCounts)
print("ALL counts")
print(length(allCounts))
transCounts <- lapply(seq(1,length(allCounts)), function(i) lapply(allCounts[[i]], function(j) j/as.numeric(bias[i])))
if (length(labels) <= 8) {
  samplePal <- brewer.pal(9,"Set1")
  tmp <- samplePal[6]
  samplePal[6] <- samplePal[8]
  samplePal[8] <- tmp
} else {
  samplePal <- rainbow(length(labels))
}


AssignColors <- function(name, categories, pal, n=1){
  for (i in seq(1,length(categories))) {
    if (length(grep(categories[i], name)) > 0) {
      print(sprintf("i: %d mod: %d,  %d\n", i, n, i%%n))
      return(i %% n)
    }
  }
}

#
# set up plotting colors
#
n <- names(transCounts[[1]])

#labels <- c("ABC9","ABC10","ABC11", "CHM1")

print("files")
print(files)
print("colors")
print(samplePal)
print(length(samplePal))
cols <- samplePal[unlist(sapply(seq(1,length(files)), function(i) AssignColors(files[i], labels, samplePal, length(samplePal))  ) )]

print(cols)

print("labels")
print(labels)
labelCols <- samplePal[unlist(sapply(seq(1,length(labels)), function(i) AssignColors(labels[i], labels, samplePal)))]
print("labelcols")
print(labelCols)
print("drawing plots")
groupIndices <- sapply(seq(1,length(files)), function(i) AssignColors(files[i], labels, samplePal))


#boxesFileName <- sprintf("%s.boxes.pdf",args$output)
#
#pdf(boxesFileName)
#sapply(seq(1,length(n)), function(i) PlotBoxes(transCounts, n[i], groupIndices, colors=cols, legLabels=labels[1:5], legPal=labelCols[1:5]))
#dev.off()
#
                                        #limTransCounts <- lapply(seq(1,length(transCounts)), function(i) lapply(transCounts[[i]], function(j) sapply(j, function(x) min(x,400))))

#
# Test plotting count


if (length(roiTable) > 0) {
  roi <- lapply(seq(1,length(n)), function(i) { idx <- which(roiTable$V1 == n[i]); if (length(roiTable) > 0) { return(c(roiTable$V4[idx[1]], roiTable$V5[idx[1]])) } else { return(c(0,0))}})
} else {
  roi <- rep(c(NA, NA), length(n))
}



coverageFileName <- sprintf("%s.coverage.pdf",args$output)

#PlotCount(transCounts, positions[[i]]-rep(positions[[i]][1],length(positions[[i]])), n[i], colors=cols, legLabels=labels, legPal=labelCols, rgn=roi[[i]])
pdf(coverageFileName)
sapply(seq(1,length(n)), function(i) PlotCount(transCounts, positions[[i]]-rep(positions[[i]][1],length(positions[[i]])), n[i], colors=cols, legLabels=labels, legPal=labelCols, rgn=roi[[i]]))
dev.off()

#coverageFileName <- "test.coverage.pdf"
#pdf(coverageFileName)
#PlotCount(limTransCounts, positions[[14]], n[14], colors=cols, legLabels=labels, legPal=labelCols)
#dev.off()
#

#log2FileName <- sprintf("%s.log2.pdf",args$output)
#pdf(log2FileName)
#sapply(seq(1,length(n)), function(i) PlotLog2(transCounts, positions[[i]], n[i], colors=cols, legLabels=labels, legPal=labelCols))
#dev.off()
#
cat(sprintf("%s\n",n),file="Titles.txt", sep='')

#sprintf("Output is in %s/%s",getwd(),boxesFileName)
sprintf("Output is in %s/%s",getwd(),coverageFileName)
#sprintf("Output is in %s/%s",getwd(),log2FileName)
