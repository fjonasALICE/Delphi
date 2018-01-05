#!/bin/bash

#SBATCH -D /gluster2/h_popp01/Delphi/
#SBATCH -o /gluster2/h_popp01/Delphi/log/py8_%j.eo.log
#SBATCH -e /gluster2/h_popp01/Delphi/log/py8_%j.eo.log

if [ "$#" -lt "3" ];
then
    echo -e "Need at least first 3 arguments of:\n[\"MB\",\"MBVeto\",\"JJ\",\"PromptPhoton\",\"WeakBoson\"]\n[number of events per pthatbin]\n[cm energy in GeV]\n[\"fullEvents\",\"noMPI\",\"noHadro\",\"noMPInoHadro\",\"noShower\"]\n[renormScaleFac]\n[factorMultFac]\n[beta_boost_z]\n[pdfA]\n[pdfB]"
    exit $?
fi

#SBATCH -J Py8
##SBATCH -p long
##SBATCH -p test
##SBATCH --exclude=node33
##SBATCH --exclude=node31
##SBATCH --exclude=node30
 
####SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000


# Something to execute
#---------------------------------------------------------------------
if [ $# == 3 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$SLURM_ARRAY_JOB_ID
elif [ $# == 4 ]; # renorm. factor
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$SLURM_ARRAY_JOB_ID
elif [ $# == 5 ]; # factor. factor
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$SLURM_ARRAY_JOB_ID
elif [ $# == 6 ]; # pythia switch for noMPI, noHadro, noMPInoHadro, noShower or fullEvents(default)
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_$SLURM_ARRAY_JOB_ID
elif [ $# == 7 ]; # with boost along z axis for asymmetric collision like pPb
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$SLURM_ARRAY_JOB_ID
elif [ $# == 8 ]; # with external PDF for beam A + B
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$8_$SLURM_ARRAY_JOB_ID
elif [ $# == 9 ]; # with external PDF for beam B only
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$8_$9_$SLURM_ARRAY_JOB_ID
fi

if [ ! -d $dirBaseName ];
then
    mkdir -p $dirBaseName
fi

if [ ! -f $dirBaseName/pythia_photons ];
then
    cp pythia_photons $dirBaseName
fi

if [ ! -f $dirBaseName/pythia_photons.cpp ];
then
    cp pythia_photons.cpp $dirBaseName
fi

if [ ! -f $dirBaseName/hendrikshelper.cxx ];
then
    cp hendrikshelper.cxx $dirBaseName
fi

sleep 1
cd $dirBaseName
sleep 1

STR1="time ./pythia_photons $SLURM_ARRAY_TASK_ID $1 $2 $3 $4 $5 $6 $7 $8 $9"

echo $STR1
eval $STR1

exit $?
