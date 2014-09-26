#!/usr/bin/env bash

rgn=`qstat -j $1 | grep rgn | awk -F'/' '{print $2}'`
#region is now a command
#cmd.rgn_2469.sh
rgndir=`echo $rgn | perl -pe 's/cmd//g'`
echo $j $rgndir

