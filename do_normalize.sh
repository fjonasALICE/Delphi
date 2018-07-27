#!/bin/bash

if [ "$#" -lt "1" ];
then
    echo "For normalization per event give path to root file as argument. Aborting"
    exit $?
fi

EXECSTRING="root -l /gluster2/h_popp01/Delphi/macros/normalize_weightSum.C'(\"${1}\")'"
echo $EXECSTRING
eval $EXECSTRING

exit
