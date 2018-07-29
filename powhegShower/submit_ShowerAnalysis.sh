#!/bin/bash

if [ "$#" -lt 7 ]; then
    echo "Give at least eight parameters:\n\
1: energy in GeV for each of the beams\n\
2: the LHAPDF ID of the PDF#1\n\
3: the LHAPDF ID of the PDF#2\n\
4: minimum bornkt in GeV\n\
5: born suppression factor in GeV\n\
6: enhanced radiation factor (0=off)\n\
7: noMPI to deactivate MPI (optional)\n\
To check which lhe files are available for showering do: sh list_available_lhef.sh\n\
exiting...\n"
    exit $?
fi

#SBATCH -J pwhgShower
#SBATCH -t 0-01:59:00
#SBATCH -p main

#Number of subjobs to start
###SBATCH --array=1-250
###SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000

#SBATCH -D /lustre/nyx/alice/users/hpoppenb/Delphi/powhegShower/

# Redirect error/output stream to a file
#SBATCH -o /lustre/nyx/alice/users/hpoppenb/Delphi/powhegShower/log/shower_%j.eo.log
#SBATCH -e /lustre/nyx/alice/users/hpoppenb/Delphi/powhegShower/log/shower_%j.eo.log

# prepare output directory
#---------------------------------------------------------------------
BEAME=${1}
PDF1=${2}
PDF2=${3}
BKTMIN=${4}
BSUP=${5}
RADFAC=${6}
NOMPI=${7}

OUTPUTDIR=beamE${BEAME}/pdf${PDF1}_${PDF2}/bktmin${BKTMIN}_bsup${BSUP}_radfac${RADFAC}
BASESTRING=beamE${BEAME}_pdf${PDF1}_${PDF2}_bktmin${BKTMIN}_bsup${BSUP}_radfac${RADFAC}
BASEDIR=$OUTPUTDIR

# change output dir if MPI deactivated
if [[ "${NOMPI}" == "noMPI"  ]]
then
    OUTPUTDIR+="_noMPI"
    BASESTRING+="_noMPI"
fi
wait


if [ ! -d ${OUTPUTDIR} ];
then
    mkdir -p ${OUTPUTDIR}
fi

if [ ! -f ${OUTPUTDIR}/shower.conf ];
then
    cp shower.conf ${OUTPUTDIR}
fi

if [ ! -f ${OUTPUTDIR}/ShowerAnalysis ];
then
    cp ShowerAnalysis ${OUTPUTDIR}
fi

if [ ! -f ${OUTPUTDIR}/ShowerAnalysis.cpp ];
then
    cp src/ShowerAnalysis.cpp ${OUTPUTDIR}
fi

wait
# execute
#---------------------------------------------------------------------
sleep 1
cd ${OUTPUTDIR}
sleep 1

# change config file for options
if [[ "${NOMPI}" == "noMPI"  ]]
then
    sed -i -e 's/#PartonLevel:MPI = off/PartonLevel:MPI = off/g' shower.conf
fi
wait

INDEX1=$((SLURM_ARRAY_TASK_ID*2-1))
INDEX2=$((SLURM_ARRAY_TASK_ID*2))

POWHEGPREFIX=/lustre/nyx/alice/users/hpoppenb/POWHEG-BOX-V2

STR1="time ./ShowerAnalysis $INDEX1.root ${POWHEGPREFIX}/directphoton/lhefGen/lhef/${BASEDIR}/${INDEX1}/pwgevents-0001.lhe &"
STR2="time ./ShowerAnalysis $INDEX2.root ${POWHEGPREFIX}/directphoton/lhefGen/lhef/${BASEDIR}/${INDEX2}/pwgevents-0001.lhe"
sleep 2
echo $STR1
eval $STR1
sleep 2
echo $STR2
eval $STR2
wait

cd -
# create a merge script for the produced root files
if [[ ! -f do_hadd_${BASESTRING}.sh ]]; then
    touch do_hadd_${BASESTRING}.sh
    cat <<EOF >do_hadd_${BASESTRING}.sh
hadd mergedOutput/${BASESTRING}_merged.root ${OUTPUTDIR}/*.root 
EOF
fi

exit $!
