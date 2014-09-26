clist=$(shell cat $(CONSENSI))
renamed=$(subst :,Z,$(clist))
fasta=$(addsuffix .fasta,$(renamed))

%.fasta:
	$(eval region=$(subst Z,:,$*))
	echo $(region)
	echo "module load python/latest" > commands/$*.csh
	echo "module load numpy/latest" >> commands/$*.csh
	echo "~mchaisso/software/bin/quiver /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio/pacbio.10x.sorted.cmp.h5 -j 4 --referenceWindows $(region) -o $@ -r /net/eichler/vol2/eee_shared/assemblies/hg19/ucsc.hg19.fasta  " >> commands/$*.csh 	
	qsub -S /bin/bash -V -cwd -sync y -pe orte 1 -l h_rss=4G commands/$*.csh 

all:$(fasta)
