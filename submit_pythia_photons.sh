#!/bin/bash

#SBATCH -J softQCDtest
#SBATCH -p long
#SBATCH -D /gluster2/h_popp01/Delphi
#SBATCH -o /gluster2/h_popp01/Delphi/log/shower_%j.eo.log
#SBATCH -e /gluster2/h_popp01/Delphi/log/shower_%j.eo.log

####SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000



# Something to execute
#---------------------------------------------------------------------

if [ ! -d pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID ];
then
    mkdir -p pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID/pythia_photons ];
then
    cp pythia_photons pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID/haddav ];
then
    cp haddav pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID
fi

sleep 1
cd pythia_photons_$1_$2_$3_$4_$5_$6_$SLURM_ARRAY_JOB_ID/
sleep 1

STR1="time ./pythia_photons $SLURM_ARRAY_TASK_ID $1 $2 $3 $4 $5 $6 "

echo $STR1
eval $STR1

exit $?
