HOMEDIR = /net/eichler/vol5/home/mchaisso

SAMTOOLS=$(HOMEDIR)/software/samtools/samtools-0.1.18/

all: findUnique guaranteeUnique queryBam queryFasta countRepeatComposition filterlc nloc pmme

BLASR_COMMON =  /net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/common

findUnique: FindUnique.cpp
	g++ -g -O3 FindUnique.cpp -I $(BLASR_COMMON) -o findUnique

countRepeatComposition: CountRepeatComposition.cpp
	g++ -g CountRepeatComposition.cpp -I $(BLASR_COMMON) -o countRepeatComposition -static

guaranteeUnique: GuaranteeUnique.cpp
	g++ -g -O3 GuaranteeUnique.cpp -I $(BLASR_COMMON) -o guaranteeUnique -lpthread


queryBam: QueryBam.cpp
	g++ -O3  QueryBam.cpp -I $(BLASR_COMMON)  -I  /net/eichler/vol5/home/mchaisso/projects/ExomeAssembly/code/simplegraph  -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -lpthread -o $@ 

queryFasta: QueryFasta.cpp
	g++ -O3 QueryFasta.cpp -I $(BLASR_COMMON)  -I  /net/eichler/vol5/home/mchaisso/projects/ExomeAssembly/code/simplegraph -lm -L$(HOMEDIR)/software/lib -lz -lpthread -o $@ -static

filterlc: FilterLowComplexity.cpp
	g++ -O3 $^ -I $(BLASR_COMMON) -I  /net/eichler/vol5/home/mchaisso/projects/ExomeAssembly/code/simplegraph -lm -L$(HOMEDIR)/software/lib -lz -lpthread -o $@ -static


nloc: PrintNLocations.cpp
	g++ -O3 $^ -I $(BLASR_COMMON) -lpthread -o $@ -static

pmme: PrintMostMaskedEntries.cpp
	g++ -O3 $^ -I $(BLASR_COMMON) -lpthread -o $@ -static
