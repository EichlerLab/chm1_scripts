#!/usr/bin/env bash

mkdir -p $3
if [ ! -s $2  ] ; then

exit 0
fi


if [ $1 = "repeatmasker" ]; then
RepeatMasker -xsmall -dir $3 -pa 8 $2
else
echo "will link  " $2 $3
cp $2 $3/
b=$(basename $2)
cd $3; $HOME/software/bin/censor.ncbi $b -lib hum -s -simple -classify ALU
fi
