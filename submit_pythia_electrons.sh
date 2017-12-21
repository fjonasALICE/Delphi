#!/bin/bash

##SBATCH -J PyElec
#SBATCH -p long
##SBATCH -p test
##SBATCH --exclude=node33
##SBATCH --exclude=node31
##SBATCH --exclude=node30
 
#SBATCH -D /gluster2/h_popp01/Delphi/
#SBATCH -o /gluster2/h_popp01/Delphi/log/pyElectron_%j.eo.log
#SBATCH -e /gluster2/h_popp01/Delphi/log/pyElectron_%j.eo.log

####SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000




# Something to execute
#---------------------------------------------------------------------

if [ ! -d pythia_electrons_$1_$SLURM_ARRAY_JOB_ID ];
then
    mkdir -p pythia_electrons_$1_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f pythia_electrons_$1_$SLURM_ARRAY_JOB_ID/pythia_electrons ];
then
    cp pythia_electrons pythia_electrons_$1_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f pythia_electrons_$1_$SLURM_ARRAY_JOB_ID/pythia_electrons.cpp ];
then
    cp pythia_electrons.cpp pythia_electrons_$1_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f pythia_electrons_$1_$SLURM_ARRAY_JOB_ID/hendrikshelper.cxx ];
then
    cp hendrikshelper.cxx pythia_electrons_$1_$SLURM_ARRAY_JOB_ID
fi

sleep 1
cd pythia_electrons_$1_$SLURM_ARRAY_JOB_ID/
sleep 1

STR1="time ./pythia_electrons $SLURM_ARRAY_TASK_ID $1"

echo $STR1
eval $STR1

exit $?
