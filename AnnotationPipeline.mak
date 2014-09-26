PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
MASKER?=repeatmasker
alignments?=alignments.sam

all: $(alignments) insertions.bed\
	deletions.bed\
	insertion/insertions.fasta\
	deletion/deletions.fasta\
	insertion/rm/insertions.fasta.out \
	deletion/rm/deletions.fasta.out \
	deletion/deletions.annotated.bed\
	insertion/insertions.annotated.bed\
	insertion/L1.bed\
	deletion/L1.bed

#
# First align the fasta to the reference.
#


$(alignments): $(INPUT_FASTA)
	-cp -f $(alignments) $(alignments).bak
	$(BLASR) $(INPUT_FASTA) /var/tmp/mchaisso/ucsc.hg19.fasta -sa /var/tmp/mchaisso/ucsc.hg19.fasta.sa -ctab /var/tmp/mchaisso/ucsc.hg19.fasta.ctab -affineAlign -affineOpen 100 -affineExtend 0 -insertion 0 -deletion 0 -nproc 12 -out $@ -bestn 1 -sam  -maxMatch 50 


#
# Next print the gapped sequences, split into insertion and deletion, and simplify redundant (duplicate) gaps.
#
gaps.bed: $(alignments)
	$(PBS)/PrintGaps.py /var/tmp/mchaisso/ucsc.hg19.fasta $(alignments) --condense 20 --tsd 20 --outFile gaps.bed

insertions.bed: gaps.bed
	grep insertion gaps.bed | bedtools sort | $(PBS)/rmdup.py /dev/stdin $@ --leftjustify --window 100
	$(PBS)/FilterEventListByAssembledTarget.py $@ $@.tmp
	/bin/mv $@.tmp $@
	grep -s -v "chrY" $@ > $@.tmp
	/bin/mv -f $@.tmp $@

deletions.bed: gaps.bed
	grep deletion gaps.bed | bedtools sort | $(PBS)/rmdup.py /dev/stdin $@ --leftjustify --window 100
	$(PBS)/FilterEventListByAssembledTarget.py $@ $@.tmp
	/bin/mv -f $@.tmp $@
	grep -s -v "chrY" $@ > $@.tmp
	/bin/mv -f $@.tmp $@


#
# First step to annotation is to create a fasta file.
#

insertion/insertions.fasta: insertions.bed
	-mkdir -p insertion
	$(PBS)/GapBedToFasta.py insertions.bed insertion/insertions.fasta

deletion/deletions.fasta: deletions.bed
	-mkdir -p deletion
	$(PBS)/GapBedToFasta.py deletions.bed deletion/deletions.fasta

#
# Next, the repeat masked files are generated.
#

insertion/rm/insertions.fasta.out:  insertion/insertions.fasta
	$(PBS)/MaskFasta.sh $(MASKER) insertion/insertions.fasta insertion/rm

deletion/rm/deletions.fasta.out:  deletion/deletions.fasta
	$(PBS)/MaskFasta.sh $(MASKER) deletion/deletions.fasta deletion/rm


deletion/deletions.annotated.bed: deletion/rm/deletions.fasta.out
	$(PBS)/AnnotateGapBed.py deletions.bed $@ deletion/rm/deletions.fasta.out deletion/rm/deletions.fasta.masked
#	$(PBS)/AnnotateGapBedWithCensorMap.py deletions.bed deletion/rm/deletions.fasta.map $@

	grep NONE deletion/deletions.annotated.bed > deletion/deletions.NONE.bed
	grep -v NONE deletion/deletions.annotated.bed > tmp
	/bin/mv tmp deletion/deletions.annotated.all.bed
	cat deletion/deletions.annotated.all.bed | awk '{ if ($$10 < 0.70) print}' > deletion/deletions.partial_masked.bed
	cat deletion/deletions.annotated.all.bed | awk '{ if ($$10 >= 0.70) print}' > $@


insertion/insertions.annotated.bed: insertion/rm/insertions.fasta.out
	$(PBS)/AnnotateGapBed.py insertions.bed $@ insertion/rm/insertions.fasta.out insertion/rm/insertions.fasta.masked
#	$(PBS)/AnnotateGapBedWithCensorMap.py insertions.bed insertion/rm/insertions.fasta.map $@
	grep NONE insertion/insertions.annotated.bed > insertion/insertions.NONE.bed
	grep -v NONE insertion/insertions.annotated.bed > tmp
	/bin/mv tmp insertion/insertions.annotated.all.bed
	cat insertion/insertions.annotated.all.bed | awk '{ if ($$10 < 0.50) print}' > insertion/insertions.partial_masked.bed
	cat insertion/insertions.annotated.all.bed | awk '{ if ($$10 >= 0.50) print}' > $@

#
# The resulting repeat masked files should have mostly sequence annotated as repeat inside of it.
#

insertion/L1.bed: insertion/insertions.annotated.bed
	$(PBS)/PrintUniqueEvents.sh insertion/insertions.annotated.bed insertion
	rm -rf insertion/[0-9]r*
	mkdir insertion/full
	mv insertion/ins* insertion/full

deletion/L1.bed: deletion/deletions.annotated.bed
	$(PBS)/PrintUniqueEvents.sh deletion/deletions.annotated.bed deletion
	rm -rf deletion/[0-9]*
	mkdir deletion/full
	mv deletion/del* deletion/full
