#!/bin/bash

PDFS="CT14nlo EPPS16nlo_CT14nlo_Pb208" # test NLO version of pdf

BOOST=0.435 # shift center of mass system

STR="sbatch -p long --array=1-10 submit_PythiaAnalysis.sh JJ 1000 8160 fullEventsMonsh 1.00 1.00 ${BOOST} ${PDFS}"

eval $STR