#!/usr/bin/env bash
if [ $# -lt 2 ]; then
  echo "usage: AssembleRegion.sh alignemnts.bam region"
  exit
fi

if [ $# -eq 3 ]; then
  region=$3
else
	region='region'
fi

if [ $# -eq 4 ]; then
		slop="--slop $4"
else
		slop=""
fi


module load python/2.7.3

~mchaisso/projects/PacBioSequencing/scripts/RegionToFasta.py $1 $2 --out reads.fasta  --max 10000 --subsample $slop
if [ $? -eq 0 ]; then
	if [ -e $region/$region.ctg.fasta ]; then
			exit 0
	fi
	rm -rf $region/[0-9]\-*
#	rm -rf $region/region*

  ~mchaisso/projects/PacBioSequencing/scripts/FastaToFakeFastq.py reads.fasta reads.fastq
  ~mchaisso/software/source/celera-svn/wgs/Linux-amd64/bin/fastqToCA  -libraryname $region -technology pacbio-raw -reads reads.fastq > reads.frg
  ~mchaisso/software/wgs-8.1/Linux-amd64/bin/runCA -p $region -d $region ovlErrorRate=0.40 utgGraphErrorRate=0.40 cnsErrorRate=0.40 cgwErrorRate=0.40 unitigger=bogart obtErrorRate=0.30  reads.frg ovlThreads=8
	cp $region/9-terminator/$region.ctg.fasta $region/
  # 
	# By default, call consenus on every assembly.
	# 
	mkdir -p /var/tmp/mchaisso
	source ~mchaisso/scripts/setup_pacbio.sh
	module load python/2.7.3
  ~mchaisso/projects/PacBioSequencing/scripts/RegionToConsensus.py $1 --region $2 --delta 10000 --tmpdir /var/tmp/mchaisso  --reference $region/$region.ctg.fasta --consensus $region/$region.ctg.consensus.fasta  --p5c3
  rm -rf $region/[0-9]\-*
	perl -pi -e "s/\|/_/g" $region/$region.ctg.consensus.fasta
else
  echo "Error selecting reads. "  $? 
fi
