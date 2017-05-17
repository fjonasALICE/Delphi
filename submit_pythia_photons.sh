#!/bin/bash

module load /cvmfs/it.gsi.de/modulefiles/compiler/gcc/6.1.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/alice.gsi.de/alicesw/root/v5-34-30/inst/lib

#SBATCH -J shower_photon
#SBATCH -t 1-01:59:00
#SBATCH -p long
##SBATCH -t 0-07:59:00
##SBATCH -p main

#Number of subjobs to start
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000

#SBATCH -D /lustre/nyx/alice/users/hpoppenb/POWHEG-BOX-V2-Rev3143/powheg-directphoton/pythia/R_gamma

# Redirect output stream to a file
#SBATCH -o  /lustre/nyx/alice/users/hpoppenb/POWHEG-BOX-V2-Rev3143/powheg-directphoton/pythia/R_gamma/log/shower_%j.eo.log
#SBATCH -e  /lustre/nyx/alice/users/hpoppenb/POWHEG-BOX-V2-Rev3143/powheg-directphoton/pythia/R_gamma/log/shower_%j.eo.log

# Function to print log messages
_log() {
  local format='+%Y/%m/%d-%H:%M:%S'
  echo [`date $format`] "$@"
}

_log Job $SLURM_JOB_ID \($SLURM_JOB_NAME\)
_log Submitted from $SLURM_SUBMIT_DIR
_log Running on $USER@`hostname`:$PWD \($SLURM_JOB_NODELIST\)


# Something to execute
#---------------------------------------------------------------------

if [ ! -d monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID ];
then
    mkdir -p monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID/monash_shower_photons ];
then
    cp monash_shower_photons monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID
fi

#if [ ! -f monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID/monash_shower_photons.conf ];
#then
#    cp monash_shower_photons.conf monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID
#fi

if [ ! -f monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID/haddav ];
then
    cp haddav monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID
fi

if [ ! -f monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID/haddav.sh ];
then
    cp haddav.sh monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID
fi

#ST1=$2
#ST2=$3

sleep 1
cd monash_shower_photons_$2_$3_$SLURM_ARRAY_JOB_ID/
sleep 1

INDEX1=$((SLURM_ARRAY_TASK_ID*2-1))
INDEX2=$((SLURM_ARRAY_TASK_ID*2))


#if [ ($((!SLURM_ARRAY_TASK_ID % 2)) eq 0 ];
#then
#    sed -ie 's/SigmaProcess:renormMultFac/SigmaProcess:renormMultFac = '${ST1}' #/g' shower_LO.conf
#    sleep 2
#    sed -ie 's/SigmaProcess:factorMultFac/SigmaProcess:factorMultFac = '${ST2}' #/g' shower_LO.conf
#    sleep 2
#fi


STR1="time ./monash_shower_photons $INDEX1.root $1 &"
sleep 1
STR2="time ./monash_shower_photons $INDEX2.root $1"

echo $STR1
eval $STR1

sleep 1

echo $STR2
eval $STR2


# The process ID of the last spawned child process
child=$!
_log $0 with PID $child
# Wait for the child to finish
wait $child
# Exit signal of the child process
state=$?


_log Finishing with $state
# Propagate last signal to the system
exit $state
