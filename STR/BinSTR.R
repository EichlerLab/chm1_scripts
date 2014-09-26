 library(getopt)

options <- matrix(c("bedfile", "b", 0, "character",
                    "fai", "f", 0, "character",
                    "size", "s",0, "integer"), byrow=T, ncol=4)

args <- getopt(options)

cc <- c("character", "integer", rep("NULL", 2), "integer")
cn <- c("chrom", "pos", rep("NULL", 2), "length")

bedTable <- read.table(args$bedfile, header=F, colClasses=cc, col.names=cn)

faiTab <- read.table(args$fai, header=F)
chrs <- faiTab$V1


nChr <- length(chrs)

StoreSTR <- function(bedfile, chr) {
  i <- which(bedfile$chrom==chr)
  return(bedfile[i,])
}

strByChr <- lapply(seq(1,nChr), function(i) StoreSTR(bedTable, chrs[i]))

countsByChr <- lapply(seq(1,nChr), function(i) hist(rep(strByChr[[i]]$pos, strByChr[[i]]$length), breaks=c(seq(0,faiTab$V2[i], by=args$size),faiTab$V2[i]), plot=F))

for (i in seq(1,nChr)) {
  binStarts <- c(seq(0,faiTab$V2[i], by=args$size), faiTab$V2[i])
  for (j in seq(1,length(countsByChr[[i]]$mids))) {
    cat(sprintf("%s\t%d\t%d\t%g\n", chrs[i], binStarts[j], binStarts[j+1], countsByChr[[i]]$counts[j]))
  }
}



#print(chrs)


