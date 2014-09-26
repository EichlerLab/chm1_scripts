input=$(wildcard */input.fofn)

base=$(dir $(input))
assemblies=$(addsuffix falcon/all_tigs.fa, $(base))

all: $(assemblies)



%/falcon/all_tigs.fa: 
	cd $* && make -f /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/falcon/RunFalconAssembly.mak 

