#!/bin/bash

#SBATCH -D /gluster2/h_popp01/Delphi/
#SBATCH -o /gluster2/h_popp01/Delphi/log/py8_%j.eo.log
#SBATCH -e /gluster2/h_popp01/Delphi/log/py8_%j.eo.log

if [ "$#" -lt "3" ];
then
    echo -e "Need at least first 3 arguments of:\n\
[\"MB\",\"MBVeto\",\"JJ\",\"GJ\",\"WeakBoson\"]\n\
[number of events per pthatbin]\n\
[cm energy in GeV]\n\
[\"fullEvents\",\"noMPI\",\"noHadro\",\"noMPInoHadro\",\"noShower\"]\n\
[renormScaleFac]\n\
[factorMultFac]\n\
[beta_boost_z]\n\
[pdfA]\n\
[pdfB]\n"
    exit $?
fi

#SBATCH -J Py8
#SBATCH -p long
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

PROCESS=$1
NEVENTS=$2
CMENERGY=$3
SHOWEROPT=$4
RENSCALE=$5
FACSCALE=$6
BOOSTZ=$7
PDF1=$8
PDF2=$9

DATE=$(date +%F)

DIRBASE=py8events_${CMENERGY}GeV

# Something to execute
#---------------------------------------------------------------------
if [ $# == 3 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 4 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 5 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_${RENSCALE}_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 6 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_${RENSCALE}_${FACSCALE}_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 7 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_${RENSCALE}_${FACSCALE}_betaZ${BOOSTZ}_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 8 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_${RENSCALE}_${FACSCALE}_betaZ${BOOSTZ}_${PDF1}_jobID$SLURM_ARRAY_JOB_ID
elif [ $# == 9 ];
then
    OUTDIR=${DIRBASE}/${DATE}_${PROCESS}_${NEVENTS}ev_${SHOWEROPT}_${RENSCALE}_${FACSCALE}_betaZ${BOOSTZ}_${PDF1}_${PDF2}_jobID$SLURM_ARRAY_JOB_ID
fi

if [ ! -d $OUTDIR ];
then
    mkdir -p $OUTDIR
fi

if [ ! -f $OUTDIR/PythiaAnalysis ];
then
    cp PythiaAnalysis $OUTDIR
fi

if [ ! -f $OUTDIR/PythiaAnalysis.cpp ];
then
    cp src/PythiaAnalysis.cpp $OUTDIR
fi

if [ ! -f $OUTDIR/PythiaAnalysisHelper.cxx ];
then
    cp src/PythiaAnalysisHelper.cxx $OUTDIR
fi

sleep 1
cd $OUTDIR
sleep 1

STR1="time ./PythiaAnalysis test.root ${PROCESS} ${NEVENTS} ${CMENERGY} ${SHOWEROPT} ${RENSCALE} ${FACSCALE} ${BOOSTZ} ${PDF1} ${PDF2}"

echo $STR1
eval $STR1

# create a merge script for the produced root files and to normalize histos per event
ROOTFILENAME="merged_${PROCESS}.root"

if [[ ! -f do_merge_normalize.sh ]]; then
    touch do_merge_normalize.sh
    cat << EOF >do_merge_normalize.sh
hadd ${ROOTFILENAME} *.root
wait
root -l -q -b -x './macros/normalize_weightSum.C("${ROOTFILENAME}")'
EOF
fi
 
exit $?
