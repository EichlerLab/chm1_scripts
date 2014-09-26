asm=$(wildcard rgn_*)
regionfile=$(addsuffix /region.txt, $(asm))
image=$(addsuffix .png, $(asm))
all: $(image)


%.png:
#	$(eval region=$(shell cat $*/region.txt))
	module load python/2.7.3 &&\
  $(PBS)/DotPlot.py --query $*/region/region.ctg.consensus.fasta --target /net/eichler/vol2/eee_shared/assemblies/hg19/ucsc.hg19.fasta --thin --regionFile  $*/region.txt  --savefig $@ --nolegend --matches dot:13
