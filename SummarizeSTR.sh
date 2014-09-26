#!/usr/bin/env bash
a=`cut -f 5 $1 | awk '{ if ($1 >= 10) print;}' | stats.py  | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "STR-all $a"
a=`cut -f 5 $1 | awk '{ if ($1 >= 10 && $1 < 50) print;}' | stats.py  | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "STR-10 $a"
a=`cut -f 5 $1 | awk '{ if ($1 >= 50) print;}' | stats.py  | tr '\n' ' ' | cut -f 6,7,9,10 -d ' '`
echo "STR-50 $a"
