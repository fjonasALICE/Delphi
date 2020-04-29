#!/bin/bash

PDFS="CT14nlo EPPS16nlo_CT14nlo_Pb208" # test NLO version of pdf
BOOST=0.435 # shift center of mass system

#sbatch --array=0-30 --job-name="testing" ./submit_PythiaAnalysis.sh JJ 1000 13000 fullEventsMonash 1.00 1.00 ${BOOST} ${PDFS}

sbatch --array=0-200 --job-name="testing" ./submit_PythiaAnalysis.sh GJ 1000 13000 fullEventsMonash 1.00 1.00 
