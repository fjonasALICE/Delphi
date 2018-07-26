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
 
##SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000

# argv[1]: process switch
# argv[2]: number of events
# argv[3]: eCM
# argv[4]: optional argument to switch off MPI, hadronization or entire shower
# argv[5]: renormMultFac
# argv[6]: factorMultFac
# argv[7]: boost in z direction (beta=v/c)
# argv[8]: external pdf beam A 
# argv[9]: external pdf beam B

# Something to execute
#---------------------------------------------------------------------
if [ $# == 3 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$SLURM_ARRAY_JOB_ID
elif [ $# == 4 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$SLURM_ARRAY_JOB_ID
elif [ $# == 5 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$SLURM_ARRAY_JOB_ID
elif [ $# == 6 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_$SLURM_ARRAY_JOB_ID
elif [ $# == 7 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$SLURM_ARRAY_JOB_ID
elif [ $# == 8 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$8_$SLURM_ARRAY_JOB_ID
elif [ $# == 9 ];
then
    dirBaseName=py8events_$3GeV/$1_$2ev_$4_$5_$6_betaZ$7_$8_$9_$SLURM_ARRAY_JOB_ID
fi

if [ ! -d $dirBaseName ];
then
    mkdir -p $dirBaseName
fi

if [ ! -f $dirBaseName/pythia_photons ];
then
    cp pythiatest $dirBaseName
fi

if [ ! -f $dirBaseName/pythia_photons.cpp ];
then
    cp src/pythia.cpp $dirBaseName
fi

if [ ! -f $dirBaseName/hendrikshelper.cxx ];
then
    cp src/hendrikshelper.cxx $dirBaseName
fi

sleep 1
cd $dirBaseName
sleep 1

STR1="time ./pythiatest $SLURM_ARRAY_TASK_ID $1 $2 $3 $4 $5 $6 $7 $8 $9"

echo $STR1
eval $STR1

exit $?
