library(getopt)


PlotElement <- function(xs,ys,l,h,maxArrow,c, d, t, xends) {
       #arrow point
  x <- c(xs+l,   xs+l-(min(l,maxArrow)), xs,   xs, xs+l-min(l,maxArrow) )
  y <- c(ys+(h/2),ys+h,                   ys+h, ys, ys)

  
  if (as.numeric(d) == 1) {
    x = (l - (x-xs)) + xs
  }
  polygon(x,y,col=c)

  # compute the y position
  #

  nc <- nchar(t)
  
  xlabpos <-  c(xs + (l/2) - ((nc*par("cxy")[1])/5), xs + (l/2) + (nc*par("cxy")[1])/5)
  if (length(xends) == 0 || length(which(xends < xlabpos[1])) == 0) {
    xends <- c(xends, xlabpos[2])
    ylevel <- length(xends)
  }
  else {
    ylevel <- min(which(xends < xlabpos[1]))
    xends[ylevel] = xlabpos[2]
  }

  text(xs+(l/2), ys-(0.5*ylevel), t, cex=0.5)
  
  return(xends)
}


require(RColorBrewer, quietly=T)

PlotRepeats <- function(y, h, maxArrow, repeats, insStart, insEnd) {

  if (length(repeats$start) == 0) {
    return()
  }
  winpal <- brewer.pal(9, "Pastel1")
  genpal <- brewer.pal(9, "Set1")

  colors=rep(genpal[9], length(repeats$start))

  genomeRep <- which(repeats$start >= insStart & repeats$start <= insEnd)
  windowRep <- setdiff(seq(1,length(repeats$start)), genomeRep)
  colors[intersect(genomeRep, grep(")", repeats$rep))] = genpal[1]
  colors[intersect(genomeRep, grep("rich", repeats$rep))] = genpal[6]
  colors[intersect(genomeRep, grep("complexity", repeats$rep))] = genpal[7]  

  colors[intersect(genomeRep, grep("L1", repeats$rep))]  = genpal[2]
  colors[intersect(genomeRep, grep("L2", repeats$rep))]  = genpal[3]
  colors[intersect(genomeRep, grep("Alu", repeats$rep))] = genpal[4]
  colors[intersect(genomeRep, grep("ERV", repeats$rep))] = genpal[5]


  colors[windowRep]  = winpal[9]
  colors[intersect(windowRep, grep(")", repeats$rep))] = winpal[1]
  colors[intersect(windowRep, grep("rich", repeats$rep))] = winpal[6]
  colors[intersect(windowRep, grep("complexity", repeats$rep))] = winpal[7]
  colors[intersect(windowRep, grep("L1", repeats$rep))]  = winpal[2]
  colors[intersect(windowRep, grep("L2", repeats$rep))]  = winpal[3]
  colors[intersect(windowRep, grep("Alu", repeats$rep))] = winpal[4]
  colors[intersect(windowRep, grep("ERV", repeats$rep))] = winpal[5]

  yLabOffsets <- c(0.5,1,1.5)
  offsetIndex = 1
  xends = c()
  for (i in seq(1,length(repeats$start))) {
    xends <-  PlotElement(repeats$start[i], y, repeats$end[i]-repeats$start[i], h, maxArrow, colors[i], repeats$dir[i], repeats$rep[i], xends)
  }

}

#
# Description of options:
#
#
# pos - position in genome of start of alignment. grab from sam

options <- matrix(c("pos",     "p", 2, "integer",
                    "insert",  "i", 2, "integer",
                    "window",  "w", 2, "integer",
                    "mei", "m", 2, "character",
                    "str", "s", 2, "character",
                    "chrom", "c", 2, "character",
                    "pt", "t", 2, "character",
                    "out", "o", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

pos <- args$pos
insert <- args$insert
window <- args$window


if (is.null(args$mei) == TRUE) {
  mei <- read.table("mei.bed", header=F)
} else {
  mei <- read.table(args$mei,header=F, col.names=c("start", "end", "dir", "rep"))
}

if (is.null(args$str) == TRUE) {
  str <- read.table("str.bed", header=F)
} else {
  str <- read.table(args$str, header=F, col.names=c("start", "end", "dir", "rep"))
}

if (is.null(args$pos) == TRUE) {
  pos <- 5447230
} else {
  pos <- as.numeric(args$pos)
}

if (is.null(args$insert) == TRUE) {
  insert <- 9223
} else {
  insert <- as.numeric(args$insert)
}

if (is.null(args$window) == TRUE) {
  window <- 2000
} else {
  window <- as.numeric(args$window)
}

if (is.null(args$pt) == TRUE) {
  pt <- NULL
} else {
  pt <- read.table(args$pt, header=F, col.names=c("start", "end", "ident"))
}




if (is.null(args$out) == TRUE) {
  x11(type="dbcairo", width=12,height=3)
} else if (length(grep("png", args$out)) == 1) {
  png(args$out, width=800,height=200,res=600,type="cairo")
} else if (length(grep("pdf", args$out)) == 1) {
  pdf(args$out, width=12,height=3)
}


par(mar=c(3,1,2,1))

maxY <- 7

if (is.null(pt) == FALSE) {
  maxY <- 9
}

plot(c(), xlim=c( 0, insert+(window*2)), ylim=c(0,maxY), axes=F, xlab="", ylab="")


genomeXStart <- c(0, window, insert+window)
genomeYStart <- rep(5,3)
genomeXEnd <- c(window, insert+window, insert+2*window)
genomeYEnd <- rep(6,3)

genomePal <- brewer.pal(3,"BuGn")
colors <- c(genomePal[1], genomePal[2], genomePal[1])
rect(genomeXStart, genomeYStart, genomeXEnd, genomeYEnd, col=colors)

mei$dir = as.numeric(mei$dir)
mei$end = as.numeric(mei$end)
mei$start = as.numeric(mei$start)
mei$rep = as.character(mei$rep)
PlotRepeats(3,1,100,mei, window, insert+window )


str$dir = as.numeric(str$dir)
str$end = as.numeric(str$end)
str$start = as.numeric(str$start)
str$rep = as.character(str$rep)

PlotRepeats(1,1,100,str, window, insert+window)
xlabPos <- c(0,window,insert+window, insert+window*2)
axis(1, at=xlabPos, label=paste(xlabPos+pos))


if (is.null(pt) == FALSE && length(pt$start) > 0) {
  text(pt$start[1], 8.3, paste("panTro4, ", pt$ident[1], "% identity"), cex=0.5, pos=4)
  sx0 <- rep(0,(length(pt$start)-1))
  sx1 <- rep(0,(length(pt$start)-1))
  for (i in seq(1,length(pt$start))) {
    rect(pt$start[i], 7, pt$end[i], 8, col=brewer.pal(9,"Blues")[ceiling(((max(pt$ident[i],.1)-.1))*9)])
    sx0[i] <- pt$end[i]
    sx1[i] <- pt$start[i+1]
  }
  segments(sx0, 7.5, sx1, 7.5)
}


if (is.null(args$chrom) == FALSE) {
  mtext(args$chrom, 1, line=2)
}

res <- dev.off()




