#!/usr/bin/env bash

SOURCE=$HOME/software/blasr_2/cpp
REFFILE=/net/eichler/vol2/eee_shared/assemblies/GRCh38/indexes/GRCh38.fasta
FOFN=$1
STUB=${FOFN%.*}
BASE=$(basename $STUB)
# Run alignments
$SOURCE/alignment/bin/blasr $FOFN $REFFILE -out /dev/stdout -sam -sa $REFFILE.sa -nproc 8 -bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -ctab $REFFILE.ctab -clipping subread | gzip > $STUB.sam.gz

mkdir -p /var/tmp/mchaisso

zcat $STUB.sam.gz | samtools view -bS - | samtools sort -o - /var/tmp/mchaisso/$BASE > $STUB.bam
samtools index $STUB.bam



