#!/bin/bash

for i in ../pythia/beamE*/*/*;
do
    ls -dl ${i}
    NFILES=`(ls -dl ${i}/* | wc -l)`
    NFILES=$((NFILES-3))
    echo "$NFILES root files\n"
done

exit $?
