files=$(wildcard rgn_*/region/region.ctg.fasta)
base=$(basename $(files))
consensus=$(addsuffix .consensus.fa,$(base))

BAM=/net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/chm1.pacbio.qual.bam

all: $(consensus)


commands: 
	-mkdir -p commands

%.consensus.fa: %.fasta commands
	$(eval regiondir=$(dir $@))
	$(eval region=$(shell cat $(regiondir)/../region.txt | tr -d "\n"))
	mkdir -p /var/tmp/mchaisso
	$(eval regionfn=$(subst :,_,$(region)))
	echo "source ~mchaisso/scripts/setup_pacbio.sh; /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/RegionToConsensus.py $(BAM) --region $(region) --reference $< --delta 3000 --consensus $@ --tmpdir /var/tmp/mchaisso " > commands/cmd.$(regionfn).sh
	qsub -e /dev/null -o /dev/null -S /bin/bash -V -cwd -sync y -pe orte 1 -l mfree=4G commands/cmd.$(regionfn).sh 

