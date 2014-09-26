all: falcon/all_tigs.fa falcon/CA/string_graph.pdf

source=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/falcon
pbsource=/net/eichler/vol5/home/mchaisso/software/pacbio/bin

falcon/all_tigs.fa: falcon/CA/pre_assembled_reads.fa
	cd falcon && . $(pbsource)/activate && $(pbsource)/falcon_overlap.py --min_len 3000 --n_core 8 --d_core 2 ../$^ > p_reads.ovp
	cd falcon && . $(pbsource)/activate && $(pbsource)/falcon_asm.py p_reads.ovp ../$^

.PHONY: input.fofn


falcon/CA/pre_assembled_reads.fa:
	mkdir -p falcon
	/bin/cp input.fofn falcon/input.fofn
	cd falcon && . $(pbsource)/activate &&  $(pbsource)/HBAR_WF.py $(source)/HBAR.cfg



falcon/CA/string_graph.pdf: falcon/all_tigs.fa
	java -cp  ~mchaisso/projects/falcon_assemblies/:/net/eichler/vol5/home/mchaisso/projects/falcon_assemblies/gephi-toolkit-demos/lib/gephi-toolkit.jar PrintGraph --graphFile falcon/string_graph.gexf --plotFile falcon/string_graph.pdf 
