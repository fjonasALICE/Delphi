#!/bin/bash

echo "HELLO"
PROCESS=$1
NEVENTS=$2
CMENERGY=$3
SHOWEROPT=$4
RENSCALE=$5
FACSCALE=$6
BOOSTZ=$7
PDF1=$8
PDF2=$9
echo "I MADE IT HERE"
time ./PythiaAnalysis test.root ${PROCESS} ${NEVENTS} ${CMENERGY} ${SHOWEROPT} ${RENSCALE} ${FACSCALE} ${BOOSTZ} ${PDF1} ${PDF2}
echo "I MADE IT HERE"

# create a merge script for the produced root files and to normalize histos per event
ROOTFILENAME="merged_${PROCESS}.root"

if [[ ! -f do_merge_normalize.sh ]]; then
    touch do_merge_normalize.sh
    cat << EOF >do_merge_normalize.sh
hadd ${ROOTFILENAME} *.root
wait
root -l -q -b -x './macros/normalize_weightSum.C("${ROOTFILENAME}")'
EOF
fi
 
exit $?
