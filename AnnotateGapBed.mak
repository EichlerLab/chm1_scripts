TARGET ?= alignments
DEST=$(addsuffix .insertion.annotated.bed, $(TARGET))
all: $(DEST)

%.insertion.annotated.bed: %.insertion.fasta %.insertion.bed %.insertion/insertion.fasta.out
	RepeatMasker -pa 8 -dir $*.insertion $*.insertion.fasta
	$(PBS)/AnnotateGapBed.py $*.insertion.bed $@ *.insertion.fasta insertion/$*.insertion.fasta.out

%.insertion/ionsertion.fasta.out: $*.insertion.fasta
	-mkdir -p $*.insertion
	RepeatMasker -pa 8 -dir $*.insertion $*.insertion.fasta
	mv $*.insertion/*.out $*.insertion/insertion.fasta.out


%.deletion/ionsertion.fasta.out: $*.deletion.fasta
	-mkdir -p $*.deletion
	RepeatMasker -pa 8 -dir $*.deletion $*.deletion.fasta
	mv $*.deletion/*.out $*.deletion/deletion.fasta.out

%.insertion.fasta: %.insertion.bed:
	$(PBS)/GapBedToFasta.py %^ $*.insertion.fasta $*.deletion.fasta


%.bed: %.fasta
	$(PBS)/PrintGaps.py /var/tmp/mchaisso/ucsc.hg19.fasta $^ --tsd 10 --condense 20 --outFile $*.all.bed
	$(PBS)/rmdup.py $*.all.bed $@
	grep "insertion" $@ > $*.insertion.bed
	grep "deletion" $@ > $*.deletion.bed


