sourceFiles = $(shell cat $(INPUT))
queryFile=$(QUERY)
distalCountTargets=$(addsuffix .counts, $(basename $(sourceFiles)))
localFiles=$(notdir $(sourceFiles))
localBase=$(basename $(localFiles))
localCountTargets=$(addsuffix .counts, $(localBase))

all: $(localCountTargets)


.PHONY: $(localFiles) $(distalCountTargets)

command_dir:
	@mkdir -p command_dir

%.counts: 
	$(eval sourceFile=$(filter %$*.gz,$(sourceFiles)))
	@echo "$(PBS)/CatSeq.sh $(sourceFile) | $(PBS)/queryFasta $(QUERY) /dev/stdin $@" > command_dir/cmd.$*.sh
	qsub -S /bin/bash -V -cwd -sync y -pe serial 4 command_dir/cmd.$*.sh
