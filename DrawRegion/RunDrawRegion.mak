bed=$(INPUT)
localbase=$(basename $(notdir $(bed)))
fasta=$(localbase).fasta
sam=$(fasta).sam
rmsk=RepeatMasked/$(fasta).out
blocks=$(fasta).block
WINDOW?=4000
hg19=/var/tmp/mchaisso/ucsc.hg19.fasta
panTro4=/net/eichler/vol2/eee_shared/assemblies/panTro4/blasrdb/panTro4.fa
all: $(fasta) $(sam) $(blocks) $(rmsk) regions.pdf

$(fasta): $(bed)
	$(PBS)/PrintInsertionContext.py $(bed) $(hg19) $(fasta) --window $(WINDOW) --min 10

$(sam): $(fasta)
	$(HOME)/software/blasr_2/cpp/alignment/bin/blasr $(fasta) $(panTro4) -sa $(panTro4).sa -ctab $(panTro4).ctab -sam -bestn 1 -out $(sam) -nproc 8 -affineAlign -affineOpen 100 -affineExtend 0 -maxMatch 100

$(blocks): $(sam)
	$(PBS)/SamToBlocks.py $(sam) $(blocks)

$(rmsk): $(fasta)
	-mkdir -p RepeatMasked
	RepeatMasker -pa 8 -dir RepeatMasked $(fasta)

regions.pdf:
	mkdir -p bedfiles
	grep ">" $(fasta) | tr -d ">" > regions.txt
	$(PBS)/DrawRegion/PrintDrawCommands.sh $(rmsk) $(blocks) $(WINDOW) < regions.txt
	pdfunite  `ls plots/*.pdf | sort -k2 -t'/' -n | tr '\n' ' '` regions.pdf



