region_dirs=$(wildcard rgn_*)
assemblies=$(addsuffix /region/region.ctg.fasta, $(region_dirs))
consensus=$(addsuffix /region/region.ctg.consensus.fasta, $(region_dirs))
bamfile?=$(ALL)/chm1.pacbio.qual.bam

all: commands $(assemblies) $(consensus)

consensus: $(consensus)


commands:
	mkdir -p commands


%region.ctg.fasta:
	mkdir -p $*
	$(eval region=$(shell cat $*../region.txt))
	$(eval cmd_name=$(subst /,_,$*))
	echo $(region)
	echo "#!/usr/bin/env bash " >  commands/cmd.$(cmd_name).sh
	echo "mkdir -p /var/tmp/mchaisso" >> commands/cmd.$(cmd_name).sh
	echo "export tmpdir=\`mktemp -d -p /var/tmp/mchaisso\`" >> commands/cmd.$(cmd_name).sh
	echo cd \$$tmpdir  >> commands/cmd.$(cmd_name).sh
	echo "module load python/2.7.3; module load numpy/latest " >> commands/cmd.$(cmd_name).sh
	echo "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/AssembleRegionFromBAMs.sh $(bamfile) $(region)" >> commands/cmd.$(cmd_name).sh
	echo mkdir -p $(PWD)/$* >> commands/cmd.$(cmd_name).sh
	echo cp \$$tmpdir/region/region.ctg.fasta $(PWD)/$@  >> commands/cmd.$(cmd_name).sh
	echo rm -rf \$$tmpdir >> commands/cmd.$(cmd_name).sh
	qsub -l h_rt=4:00:00 -e /dev/null -o /dev/null -S /bin/bash -V -cwd -sync y -pe orte 1 -l mfree=8G commands/cmd.$(cmd_name).sh


%region.ctg.consensus.fasta: %region.ctg.fasta
	mkdir -p $*
	$(eval region=$(shell cat $*../region.txt))
	$(eval base=$(basename $(basename $(notdir $@))))
	$(eval cmd_name=$(subst /,_,$*))
	echo "mkdir -p /var/tmp/mchaisso" > commands/cmd.$(cmd_name).cons.sh
	echo "export tmpdir=\`mktemp -d -p /var/tmp/mchaisso\`" >> commands/cmd.$(cmd_name).cons.sh
	echo "source ~/scripts/setup_pacbio.sh; module load python/2.7.3; module load numpy/latest; export PYTHONPATH=$(PBS):\$$PYTHONPATH; cd \$$tmpdir; /net/eichler/vol19/projects/CHM1_project/nobackups/jlhudd/gaps/bin/RegionToConsensusBAMs.py  $(bamfile) --region $(region) --consensus \$$tmpdir/region.ctg.consensus.fasta --reference $(PWD)/$^ --delta 10000; perl -pi -e \"s/\|/_/g\" $(base).consensus.fasta" >> commands/cmd.$(cmd_name).cons.sh
	echo cp \$$tmpdir/region.ctg.consensus.fasta $(PWD)/$@  >> commands/cmd.$(cmd_name).cons.sh
#	echo rm -rf \$$tmpdir >> commands/cmd.$(cmd_name).cons.sh 
	qsub -l h_rt=4:00:00 -e /dev/null -o /dev/null -S /bin/bash -V -cwd -sync y -pe orte 1 -l mfree=8G commands/cmd.$(cmd_name).cons.sh

all: $(assemblies)

consensus: $(consensus)
