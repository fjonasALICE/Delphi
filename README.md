This piece of software aims for an easier usage of Pythia8 and the Parton Showering of Powheg events.
Funtions that are often used, such as "find shower photon in event and fill to histogram", are defined in src/PythiaAnalysisHelper.h.
They are used for pythia events generated in the main directory ./ and also for powheg+shower events in ./powhegShower/  .
The shell scripts "run*sh", "submit*sh" can be used to generate events, where the pythia configuration can be given as arguments, e.g.

PDFS="NNPDF23_nlo_as_0118_qed" # test NLO version of pdf
BOOST=0.435 # shift center of mass system
STR="sbatch -p main -t 0-07:59:59 --array=1-300 submit_PythiaAnalysis.sh JJ 100000 8160 noMPI 1.00 1.00 ${BOOST} ${PDFS}"
eval $STR
