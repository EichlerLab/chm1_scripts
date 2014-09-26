require(getopt, quietly=T)
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

getBase <- function(i) {
  b=unlist(strsplit(i,"[.]"))
  return(b[1])
}

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
#  points(positions[indices], movingAverage[indices], col=color)
}

MAvgPoints <- function(counts, window) {
  windowSize <- min(window, length(counts))
  if (length(counts) == 0) {
    return(0);
  }
  movingAverage <- mav(counts, windowSize)
  indices <- which(is.na(movingAverage) == T)
  movingAverage[indices] = 0
  return(movingAverage)
}

RunMAvg <- function(data, n,window) {
  avgData <- data;
  nSeq = length(avgData)
  for (i in seq(1,nSeq)) {
    for (j in seq(1,length(n))) {
      avgData[[i]][[n[j]]] <- MAvgPoints(data[[i]][[n[j]]],window);
    }
  }
  return(avgData)
}

MAvgCounts <- function(counts, positions, name, window=10) {
  nSeq <- length(counts)
  print(sprintf("averaging %s",name))
  return(lapply(seq(1,nSeq), function(i)  MAvgPoints(counts[[i]][[name]], window)));
}


NormalizeCounts <- function(counts, name) {
  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else{  DrawCountPoints(counts[[i]][[name]], positions, colors[i], window); }})
  

}

WindowMin <- function(vec, window) {
  minVals <- sapply(seq(1,length(vec)), function(i) min(vec[i:i+window]))
  minVals[which(is.na(minVals))] <- 0
  return(minVals)
}

GetMax <- function(counts, positions, name, region, window=10) {
  nSeq <- length(counts)
  mat <- as.matrix(sapply(seq(1,nSeq), function(idx1) counts[[idx1]][[name]][region[1]:region[2]]), ncol=nSeq)
  if (dim(mat)[2] == nSeq) {
    minVals <- sapply(seq(1,nSeq), function(idx) WindowMin(mat[,idx],window))
    return(max(apply(minVals,1,max) - apply(minVals,1,min)))
  }
  else {
    return(0)
  }
}


GetMaxI <- function(counts, positions, name, window=10) {

  nSeq <- length(counts)
  mat <- as.matrix(sapply(seq(1,nSeq), function(i) transCounts[[i]][[name]]), ncol=nSeq)
  if (dim(mat)[2] == nSeq) {
    return(sapply(seq(1,nSeq), function(i) WindowMaxMin(mat[,i],20)))
  }
  else {
    return(rep(0,nSeq))
  }
}



PlotLog2 <- function(counts, positions, name, window=100, colors=c("black"), legLabels=c(""), legPal=c("black")) {
  print(sprintf("log plotting %s\n",name))
  nCounts <- length(counts)
  means <- sapply(seq(1,nCounts), function(i) mean(counts[[i]][[name]]))
  medianI <- order(means)[floor(length(means)/2)]
#  i <- 3
#  counts[[i]][[name]]/counts[[medianI]][[name]]
  
  logVals <-  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else { vals <-log(counts[[i]][[name]]/counts[[medianI]][[name]],2);  vals[which(is.na(vals))] = 0;  vals[which(is.infinite(vals))] = 0; return(vals)  }})
  
  plot(c(), ylim=c(min(sapply(logVals, min)), max(sapply(logVals, max))), xlim=c(0, max(positions)), xlab="Position in insert", ylab="Read count",main=name)
  lapply(seq(1,nCounts), function(i) {if (length(counts[[i]][[name]]) == 0) { return(0)} else { vals <-log(counts[[i]][[name]]/counts[[medianI]][[name]],2);  vals[which(is.na(vals))] = 1; DrawCountPoints(vals , positions, colors[i], window); }})

  legend("topleft", legend=legLabels, col=legPal, lty=1)

}

AssignColors <- function(name, categories, pal){
  for (i in seq(1,length(categories))) {
    print("to grep")
    print(categories[i])
    print(name)
    if (length(grep(categories[i], name)) > 0) {
      return(i)
    }
  }
}


FindRangeIndices <- function(coordinates, start, end) {
  if (length(coordinates) < 2) {
    return(c(0,0));
  }
  if (start <= coordinates[1]) {
    si <- 1
  } else {
    if (length(coordinates) > 0) {
      si <- max(which(coordinates < start))
      if (is.na(si)) {
        si <- 1
      }
    } else {
      si <- 1
    }
  }

  if (end > max(coordinates)) {
    ei <- length(coordinates)
  } else {
    t <- min(which(coordinates > end))
    ei <- max(t-1,1)
  }
  if (is.na(ei)) {
    ei <- 1
  }
  return(c(si,ei))
}

CountsRangeMean <- function(c, samp, i, s, e) {
  return(mean(c[[samp]][[i]][s:e]))
}


CountSampleMeans <- function(counts, coords, roi, i) {
  r <- FindRangeIndices(coords[[i]] - coords[[i]][1], roi[[i]][1], roi[[i]][2])
  nSamples <- length(counts)
  print(r)
  return(unlist(lapply(seq(1,nSamples), function(s) mean(counts[[s]][[i]][r[1]:r[2]]))))
}

options <- matrix(c("input", "i", 2, "character",
                    "outputbase", "o", 2, "character",
                    "norm", "n", 1, "character",
                    "nrows", "r", 1, "integer",
                    "labels", "l", 1, "character"
                    ), byrow=T, ncol=4)

args <- getopt(options)







setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/Genotyping/STR")




# testing code
#setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/STR/by_position")
#setwd("/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/analysis/CombinedQuiveredRegions/insertion/Complex")


if (is.null(args$nrows)) {
  nRows = -1
} else {
  nRows = args$nrows
}

files <- read.table(args$input)
files <- paste(files$V1)

if (is.null(args$labels)) {
  labels <- c("WEA", "EA", "SA", "AFR", "ADM")
} else {
  lt <- read.table(args$labels)
  labels <- paste(lt$V1)
  print(labels)
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

bias <- sapply(fileBase, function(f) biasTable$V2[grep(f, biasTable$V1)])
bias = as.numeric(bias)
# flatline anything that is not matched.
bias[which(is.na(bias))] = 1


nCounts <- length(allCounts)
transCounts <- lapply(seq(1,length(allCounts)), function(i) lapply(allCounts[[i]], function(j) j/as.numeric(bias[i])))
samplePal <- brewer.pal(9,"Set1")
tmp <- samplePal[6]
samplePal[6] <- samplePal[8]
samplePal[8] <- tmp



#
# set up plotting colors
#
n <- names(transCounts[[1]])

#labels <- c("ABC9","ABC10","ABC11", "CHM1")
print("files")
print(files)
cols <- samplePal[unlist(sapply(seq(1,length(files)), function(i) AssignColors(files[i], labels, samplePal)))]
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
limTransCounts <- lapply(seq(1,length(transCounts)), function(i) lapply(transCounts[[i]], function(j) sapply(j, function(x) min(x,400))))

#
# Test plotting count


coverageFileName <- sprintf("%s.coverage.pdf",args$output)
pdf(coverageFileName)
sapply(seq(1,length(n)), function(i) PlotCount(limTransCounts, positions[[i]], n[i], colors=cols, legLabels=labels, legPal=labelCols))
dev.off()


log2FileName <- sprintf("%s.log2.pdf",args$output)
pdf(log2FileName)
sapply(seq(1,length(n)), function(i) PlotLog2(transCounts, positions[[i]], n[i], colors=cols, legLabels=labels, legPal=labelCols))
dev.off()



sprintf("Output is in %s/%s",getwd(),boxesFileName)
sprintf("Output is in %s/%s",getwd(),coverageFileName)
sprintf("Output is in %s/%s",getwd(),log2FileName)


