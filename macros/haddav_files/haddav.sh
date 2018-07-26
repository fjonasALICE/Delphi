#!/bin/bash

STR="./haddav"

DIV=$((${2}/${1}))
echo "make one batch from $DIV files"

for j in `seq 1 $1`;
do
    INT1=$(( ${DIV}*(${j}-1)+1 ))
    echo $INT1
    INT2=$(( ${DIV}*${j} ))
    echo $INT2

    STR="$STR out_$j.root"
    for i in `seq $INT1 $INT2`;
    do
        STR="$STR $i.root "
    done
    STR="$STR &"
    echo $STR
    eval $STR

    sleep 2
    
    unset STR
    STR="./haddav"
    
done

exit $?
