require(RColorBrewer, quietly=T)

AddIntronLines <- function(t, exonLow,exonHigh) {
  yStart = (exonLow+exonHigh)/2
  yEnd   = exonHigh - ((exonHigh-exonLow)/4)


  for (i in seq(1, dim(t)[1]-1)) {
    midpoint=(t$V3[i]+t$V2[i+1])/2
    x0 <- c(t$V3[i], midpoint)
    y0 <- c(yStart, yEnd)
    x1 <- c(midpoint, t$V2[i+1])
    y1 <- c(yEnd, yStart)
    segments(x0, y0, x1, y1)
  }
}

DrawRepeat <- function(p, l, trackLow, trackHigh, unitSize, fillPal, borderPal) {
  #
  # draw a triangle with peak at track high, and bottom at track low
  # and with a base according a log-scale
  #

  baseWidth = l*unitSize
  halfWidth = baseWidth/2
  xpoints=c(p-halfWidth,p,p+halfWidth)
  ypoints=c(trackLow, trackHigh, trackLow)
  polygon(xpoints, ypoints, col=fillPal, border=borderPal,lwd=0.5)#,border=color)
}

DrawRepeats <- function(t, trackLow, trackHigh, unitSize, fillPal, borderPal, labelRepeat=FALSE) {
  xends <- c()
  for (i in seq(1,dim(t)[1])) {

    xs <- (t$V2[i]+t$V3[i])/2
    DrawRepeat((t$V2[i]+t$V3[i])/2, t$V3[i]-t$V2[i], trackLow, trackHigh, unitSize, fillPal[i], borderPal[i])
    if (labelRepeat == TRUE) {
      print(t)
      lab=as.character(t$V7[i])
      print(sprintf("splitting '%s'", lab))
      labTab <- table(strsplit(lab,","))
      label <- names(labTab)[which(labTab == max(labTab))[1]]

      nc <- nchar(label)
      l <- t$V3[i]-t$V2[i]
      xlabpos <-  c(xs + (l/2) - ((nc*par("cxy")[1])/5), xs + (l/2) + (nc*par("cxy")[1])/5)
      if (length(xends) == 0 || length(which(xends < xlabpos[1])) == 0) {
        xends <- c(xends, xlabpos[2])
        ylevel <- length(xends)
      }
      else {
        ylevel <- min(which(xends < xlabpos[1]))
        xends[ylevel] = xlabpos[2]
      }
      
      text(t$V2[i], (ylevel-1)*0.25+0.1, labels=sprintf("%s, %d",label, t$V3[i]-t$V2[i]))

    }
  }
}

PlotGenePlusInsertions <- function(exons, strs, expansions, xlabel, name, pal, expansion=1) {
  padding <- 0.1
  # Expansion STR track
  expLow <- .8
  expHeight <- .75
  # STR reference track
  strLow <- expLow + expHeight + padding
  strHeight <- 0.5
  # Exon track
  exonLow <- strLow + strHeight + padding
  exonHeight <- 0.75
  #
  # palette for 5'utr exon, 3'utr
  #


  par(mar=c(5,7,3,3))
  extra=max((exons$V3-exons$V2)*0.25,10000)
  xl <- c(min(exons$V2)-extra, max(exons$V3)+extra)
  plot(0,xlim=xl, ylim=c(0,exonLow+exonHeight+padding), axes=F, ylab="", xlab=xlabel, main=title)
  axis(1)

  minExonPlotWidth= ceiling((xl[2] - xl[1])/500)
  lasti <- length(pal)

  for ( i in seq(1,length(exons$V2)) ) {
    if (exons$V5[i] == "3p_utr") {
      coli = lasti-4
    }  else if (exons$V5[i] == "5p_utr") {
      coli = lasti-2
    } else {
      coli = lasti
    }

    rect(exons$V2[i], exonLow, max(exons$V3[i], exons$V2[i]+minExonPlotWidth), exonLow+exonHeight, col=pal[coli], border=pal[4])#, border=F)
  }

  AddIntronLines(exons, exonLow,exonLow+exonHeight)


  #
  # red ish, blue ish, green ish
  #
  colorKey=c("#e41a1c","#377eb8","#4daf4a")
  print(length(strs))
  if (length(strs) > 0) {
    fillPal <- rep(brewer.pal(3,"Greens")[2], length(strs$V2))
    borderPal <- rep(brewer.pal(3,"Greens")[3], length(strs$V2))
    DrawRepeats(strs, strLow, strLow+strHeight, expansion, fillPal, borderPal, length(strs$V2))
  }

  axis(2, at=c(exonLow+0.25,strLow+0.25,expLow+0.25), labels=c("Exons", "STR", "STR expansion"), tick=F, las=2)

  #
  # determine the color of each expansion
  #  exon - red
  #  intron - blue
  #

  exonPal <- brewer.pal(3,"Reds")
  intronPal <- brewer.pal(3,"Purples")
  if (length(expansions) > 0) {
    fillColors = rep(intronPal[2], length(expansions$V2))
    borColors = rep(intronPal[3], length(expansions$V2))
    for (i in seq(1,length(exons$V2))) {
      for (j in seq(1,length(expansions$V2))) {
        if ( exons$V2[i] <= expansions$V2[j] && exons$V3[i] >= expansions$V2[j]) {
          fillColors[j] = exonPal[2]
          borColors[j] = exonPal[3]
        }
      }
    }
    
    DrawRepeats(expansions, expLow, expLow+expHeight, expansion, fillColors, borColors, labelRepeat=T)
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  print("Usage: Rscript PlotGenePlusInsertions.R exonBedFileName ReferenceSTR.bed InsertedSTR.bed title outfile expansion")
  print("expansion is the ratio that the base of every triangle is expanded in the plot. For example, 10 is 10x the actual")
  print("length.")
  stop()
}


exonBedFileName = args[1]
refSTRBedFileName = args[2]
insertSTRBedFileName = args[3]
title=args[4]
outputFileName = args[5]
expansion=as.integer(args[6])


exons <- read.table(exonBedFileName, header=F)


i <- file.info(refSTRBedFileName)
if (i$size > 0) {
  refSTR <- read.table(refSTRBedFileName, header=F)
} else {
  refSTR <- c()
}


i <- file.info(insertSTRBedFileName)
if (i$size > 0) {
  insertSTR <- read.table(insertSTRBedFileName, header=F)
} else {
  insertSTR= c()
}
waitAtEnd = F
if (length(grep("png", outputFileName)) == 1) {
  png(outputFileName,width=10,height=4, type="cairo", res=600, units="in", family=X11Fonts()$Helvetica)
} else if (length(grep("pdf", outputFileName)) == 1) {
  print("writing pdf")
  pdf(outputFileName, width=10,height=4)
} else {
  waitAtEnd = T
  x11(width=10,height=4,type="dbcairo")
}

PlotGenePlusInsertions(exons, refSTR, insertSTR, exons$V1[1], title, brewer.pal(7, "Blues"), expansion)
if (waitAtEnd) {
  x <- readLines(file("stdin"),1)
}
d <- dev.off()

