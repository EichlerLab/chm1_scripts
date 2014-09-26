#!/usr/bin/env bash
if [[ "$#" -lt "2" ]]
then
  echo "usage: $0 list_of_bam_paths region [region_name]"
  exit 1
fi

if [ $# -eq 3 ]; then
    region=$3
else
    region='region'
fi

TMP_DIR=/var/tmp/`whoami`

module load python/2.7.2

# Get all reads from multiple BAMs corresponding to the given region.
INPUTS=`awk '{ print "-in "$1 }' $1 | tr '\n' ' '`
BAMTOOLS_REGION=${2/-/..}
MIN_QUALITY=20
/net/eichler/vol4/home/jlhudd/src/bamtools-2.3.0/bin/bamtools filter ${INPUTS} -region ${BAMTOOLS_REGION} -mapQuality ">=${MIN_QUALITY}" \
    | samtools view - \
    | awk '{ print ">"$1; print $10 }' > reads.fasta

if [[ "$?" -eq "0" ]]
then
    if [ -e $region/$region.ctg.fasta ]; then
        exit 0
    fi
    rm -rf $region/[0-9]\-*

    ~mchaisso/projects/PacBioSequencing/scripts/FastaToFakeFastq.py reads.fasta reads.fastq
    ~mchaisso/software/source/celera-svn/wgs/Linux-amd64/bin/fastqToCA  -libraryname $region -technology pacbio-raw -reads reads.fastq > reads.frg
    ~mchaisso/software/wgs-8.2beta/Linux-amd64/bin/runCA  -p $region -d $region ovlErrorRate=0.40 utgGraphErrorRate=0.40 cnsErrorRate=0.40 cgwErrorRate=0.40 unitigger=bogart obtErrorRate=0.30  reads.frg ovlThreads=8
    cp $region/9-terminator/$region.ctg.fasta $region/

    #
    # By default, call consenus on every assembly.
    #
    mkdir -p ${TMP_DIR}
    source ~mchaisso/scripts/setup_pacbio.sh
    export PYTHONPATH=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts:$PYTHONPATH
    python /net/eichler/vol19/projects/CHM1_project/nobackups/jlhudd/gaps/bin/RegionToConsensusBAMs.py $1 --region $2 --delta 30000 --tmpdir ${TMP_DIR}  --reference $region/$region.ctg.fasta --consensus $region/$region.ctg.consensus.fasta
    rm -rf $region/[0-9]\-*
else
    echo "Error selecting reads. "  $?
fi
