samfiles=$(wildcard m*.sam)
bedfiles=$(addsuffix .bed, $(samfiles))
%.sam.bed: %.sam
	module load python/latest; /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/DetectInversions.py --sam $< --out $@ --tests $@.fasta

all: $(bedfiles)
