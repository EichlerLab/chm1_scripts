#!/usr/bin/env bash

chr=$1
start=$2
end=$3
refdir=$4

source /net/eichler/vol5/home/mchaisso/scripts/setup_pacbio.sh
mkdir -p /var/tmp/mchaisso
$PBS/RegionToConsensus.py /net/eichler/vol20/projects/pacbio/nobackups/CHM1_from_pacbio.all/chm1.pacbio.qual.bam --region $chr:$start-$end --reference $refdir/region.ctg.fasta --consensus $refdir/reads.ctg.fasta.consensus.fa --delta 5000 --tmpdir /var/tmp/mchaisso 

