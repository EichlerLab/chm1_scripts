#!/usr/bin/env bash

/net/eichler/vol5/home/mchaisso/software/blasr_1/cpp/alignment/bin/blasr $1 $1 -bestn 30 -extend -maxExtendDropoff 100 -maxMatch 25 -out $1.sam -sam -preserveReadTitle
