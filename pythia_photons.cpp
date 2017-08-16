#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "hendrikshelper.h"

using std::cout;
using namespace Pythia8;

int main(int, char **);
int main(int argc, char **argv) {



  //--- read commandline args ----------------------------------------
  if (argc < 7) {
    printf("Invalid number of arguments\nUsage: %s [output root file] [\"dir\" of \"frag\" or \"decay\"] [number of events] [cm energy in GeV] [renormScaleFac] [factorMultFac] \nOne of the following optional is also allowed:\n You can switch off the entire parton shower with additional argument \"noShower\".\n You can switch off hadronisation with additional argument \"noHadro\".\n You can switch off MPI with additional argument \"noMPI\".\n You can switch off hadronisation and MPI with additional argument \"noMPInoHadro\".", argv[0]);
    exit(EXIT_FAILURE);
  }


  //--- define output root file --------------------------------------
  char rootFileName[1024];
  if (argc == 7)
    snprintf( rootFileName, sizeof(rootFileName), "pythia_%s_%s_%ldGeV_RenS%.2f_FacS%.2f.root", argv[1], argv[2],
              strtol(argv[4], NULL, 10), strtof(argv[5], NULL), strtof(argv[6], NULL) );
  if (argc == 8)
    snprintf( rootFileName, sizeof(rootFileName), "pythia_%s_%s_%ldGeV_RenS%.2f_FacS%.2f_%s.root", argv[1], argv[2],
              strtol(argv[4], NULL, 10), strtof(argv[5], NULL), strtof(argv[6], NULL), argv[7]);
  printf("------------------------------------------------------------------------------------------\n");
  printf("The result will be written into %s\n", rootFileName);
  printf("------------------------------------------------------------------------------------------\n");

  int nEvent = strtol(argv[3], NULL, 10); // number of events

  Pythia p;
  HendriksHelper pyHelp;

  pyHelp.Pass_Parameters_To_Pythia(p, argc, argv); // which energy, scales, optional master switches

  string pdfA = "LHAPDF6:nCTEQ15_1_1.LHgrid";
  //  string pdfB = "LHAPDF6:nCTEQ15_1_1.LHgrid";
  string pdfB = "LHAPDF6:nCTEQ15FullNuc_208_82.LHgrid";
  const bool yesBoost = true;
  p.readString("PDF:pSet = " + pdfA);
  p.readString("PDF:pSetB = " + pdfB);

  p.readString("Next:NumberCount = 100000");
  pyHelp.Set_Pythia_Randomseed(p);



  //--- Histograms ---------------------------------------------------
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const double ptMin = 0., ptMax = 100.;
  const int ptBins = 100;
  const double etaTPC = 0.9,
    etaEMCal = 0.27,
    etaPHOS = 0.12;

  // pTHatBins for direct photons
  // const int pTHatBins = 4;
  // double pTHatBin[pTHatBins+1] = {15., 20., 26., 33., 1000.};
 
  //   pTHatBins for shower photons with softQCD limit 18 GeV
  // const int pTHatBins = 13;
  // double pTHatBin[pTHatBins+1] = { 0. , 18., 21., 24., 27.,
  // 				   30., 34., 38., 43., 48.,
  // 				   56., 67., 80.,
  // 				   1000. };

  // pTHatBins for shower photons with softQCD limit 15 Gev
  const int pTHatBins = 14;
  double pTHatBin[pTHatBins+1] = { 0., 15., 18., 21.,
   				   24., 27., 30., 34., 38.,
   				   43., 48., 56., 67., 80.,
   				   1000. };

  // pTHatBins for shower photon test at lowest pt
  // const int pTHatBins = 1;
  // double pTHatBin[pTHatBins+1] = {0., 15.};
 
  printf("-----------------------\nusing %d pTHat bins\n-----------------------\n", pTHatBins);

  // controls plot to see how much is cut away with sQCD fluctuation cut
  vector <TH1D> vec_fluctCut;
    vec_fluctCut.push_back( TH1D("h_fluctCut00", "h_fluctCut00", ptBins, ptMin, ptMax) );
  for( int i = 10; i <= 30; i++ ){ 
    vec_fluctCut.push_back( TH1D(Form("h_fluctCut%02d",i), Form("h_fluctCut%02d",i), ptBins, ptMin, ptMax) );
  }

  // primary pions
  TH1D *h_pi0_yTPC = new TH1D("h_pi0_yTPC","pi0_yTPC", ptBins, ptMin, ptMax);
  TH1D *h_pi0_yEMCal = new TH1D("h_pi0_yEMCal","pi0_yEMCal", ptBins, ptMin, ptMax);
  TH1D *h_pi0_yPHOS = new TH1D("h_pi0_yPHOS","pi0_yPHOS", ptBins, ptMin, ptMax);

  // all etas without secondary correction
  TH1D *h_eta_yTPC = new TH1D("h_eta_yTPC","eta_yTPC", ptBins, ptMin, ptMax);
  TH1D *h_eta_yEMCal = new TH1D("h_eta_yEMCal","eta_yEMCal", ptBins, ptMin, ptMax);
  TH1D *h_eta_yPHOS = new TH1D("h_eta_yPHOS","eta_yPHOS", ptBins, ptMin, ptMax);

  // direct photons (consider only non-decay photons)
  TH1D *h_non_decay_photons_etaTPC = new TH1D("h_non_decay_photons_etaTPC","non_decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_non_decay_photons_etaEMCal = new TH1D("h_non_decay_photons_etaEMCal","non_decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_non_decay_photons_etaPHOS = new TH1D("h_non_decay_photons_etaPHOS","non_decay photons_etaPHOS", ptBins, ptMin, ptMax);

  // isolated photons (considers only non-decay photons)
  TH1D *h_iso_charged2GeV_R03_photons_etaTPC = new TH1D("h_iso_charged2GeV_R03_photons_etaTPC","non_decay iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R03_photons_etaEMCal","non_decay iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R03_photons_etaPHOS","non_decay iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R04_photons_etaTPC = new TH1D("h_iso_charged2GeV_R04_photons_etaTPC","non_decay iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R04_photons_etaEMCal","non_decay iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R04_photons_etaPHOS","non_decay iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R05_photons_etaTPC = new TH1D("h_iso_charged2GeV_R05_photons_etaTPC","non_decay iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R05_photons_etaEMCal","non_decay iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R05_photons_etaPHOS","non_decay iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R03_photons_etaTPC = new TH1D("h_iso_charged3GeV_R03_photons_etaTPC","non_decay iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R03_photons_etaEMCal","non_decay iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R03_photons_etaPHOS","non_decay iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R04_photons_etaTPC = new TH1D("h_iso_charged3GeV_R04_photons_etaTPC","non_decay iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R04_photons_etaEMCal","non_decay iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R04_photons_etaPHOS","non_decay iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R05_photons_etaTPC = new TH1D("h_iso_charged3GeV_R05_photons_etaTPC","non_decay iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R05_photons_etaEMCal","non_decay iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R05_photons_etaPHOS","non_decay iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R03_photons_etaTPC = new TH1D("h_iso_full2GeV_R03_photons_etaTPC","non_decay iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaEMCal = new TH1D("h_iso_full2GeV_R03_photons_etaEMCal","non_decay iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaPHOS = new TH1D("h_iso_full2GeV_R03_photons_etaPHOS","non_decay iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R04_photons_etaTPC = new TH1D("h_iso_full2GeV_R04_photons_etaTPC","non_decay iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaEMCal = new TH1D("h_iso_full2GeV_R04_photons_etaEMCal","non_decay iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaPHOS = new TH1D("h_iso_full2GeV_R04_photons_etaPHOS","non_decay iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R05_photons_etaTPC = new TH1D("h_iso_full2GeV_R05_photons_etaTPC","non_decay iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaEMCal = new TH1D("h_iso_full2GeV_R05_photons_etaEMCal","non_decay iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaPHOS = new TH1D("h_iso_full2GeV_R05_photons_etaPHOS","non_decay iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R03_photons_etaTPC = new TH1D("h_iso_full3GeV_R03_photons_etaTPC","non_decay iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaEMCal = new TH1D("h_iso_full3GeV_R03_photons_etaEMCal","non_decay iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaPHOS = new TH1D("h_iso_full3GeV_R03_photons_etaPHOS","non_decay iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R04_photons_etaTPC = new TH1D("h_iso_full3GeV_R04_photons_etaTPC","non_decay iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaEMCal = new TH1D("h_iso_full3GeV_R04_photons_etaEMCal","non_decay iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaPHOS = new TH1D("h_iso_full3GeV_R04_photons_etaPHOS","non_decay iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R05_photons_etaTPC = new TH1D("h_iso_full3GeV_R05_photons_etaTPC","non_decay iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaEMCal = new TH1D("h_iso_full3GeV_R05_photons_etaEMCal","non_decay iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaPHOS = new TH1D("h_iso_full3GeV_R05_photons_etaPHOS","non_decay iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  // decay photons
  TH1D *h_decay_photons_etaTPC = new TH1D("h_decay_photons_etaTPC","decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaEMCal = new TH1D("h_decay_photons_etaEMCal","decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaPHOS = new TH1D("h_decay_photons_etaPHOS","decay photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", 1000, ptMin, ptMax);

  TH1D *h_weightSum = new TH1D("h_weightSum","sum of weights", 1, 0., 1.);

  // organise pTHat wise histograms in vectors
  vector <TH1D*> vec_weightSum_bin;

  vector <TH1D*> vec_non_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_non_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_non_decay_photons_etaPHOS_bin;

  vector <TH1D*> vec_pi0_yTPC_bin;
  vector <TH1D*> vec_pi0_yEMCal_bin;
  vector <TH1D*> vec_pi0_yPHOS_bin;

  vector <TH1D*> vec_eta_yTPC_bin;
  vector <TH1D*> vec_eta_yEMCal_bin;
  vector <TH1D*> vec_eta_yPHOS_bin;

  vector <TH1D*> vec_iso_charged2GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_charged2GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_charged2GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_iso_charged2GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_charged2GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_charged2GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_iso_charged2GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_charged2GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_charged2GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_iso_charged3GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_charged3GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_charged3GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_iso_charged3GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_charged3GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_charged3GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_iso_charged3GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_charged3GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_charged3GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_iso_full2GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_full2GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_full2GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_iso_full2GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_full2GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_full2GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_iso_full2GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_full2GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_full2GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_iso_full3GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_full3GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_iso_full3GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_iso_full3GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_full3GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_iso_full3GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_iso_full3GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_full3GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_iso_full3GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_decay_photons_etaPHOS_bin;

  vector <TH1D*> vec_pTHat_bin;


  for(int i = 0; i < pTHatBins; i++){

    // pi0 histos
    vec_pi0_yTPC_bin.push_back( (TH1D*)h_pi0_yTPC->Clone(Form("h_pi0_yTPC_bin_%02d",i)) );
    vec_pi0_yEMCal_bin.push_back( (TH1D*)h_pi0_yEMCal->Clone(Form("h_pi0_yEMCal_bin_%02d",i)) );
    vec_pi0_yPHOS_bin.push_back( (TH1D*)h_pi0_yPHOS->Clone(Form("h_pi0_yPHOS_bin_%02d",i)) );

    // eta histos
    vec_eta_yTPC_bin.push_back( (TH1D*)h_eta_yTPC->Clone(Form("h_eta_yTPC_bin_%02d",i)) );
    vec_eta_yEMCal_bin.push_back( (TH1D*)h_eta_yEMCal->Clone(Form("h_eta_yEMCal_bin_%02d",i)) );
    vec_eta_yPHOS_bin.push_back( (TH1D*)h_eta_yPHOS->Clone(Form("h_eta_yPHOS_bin_%02d",i)) );

    // non-decay photon histos
    vec_weightSum_bin.push_back( (TH1D*)h_weightSum->Clone(Form( "h_weightSum_bin_%02d", i )) );
    vec_non_decay_photons_etaTPC_bin.push_back( (TH1D*)h_non_decay_photons_etaTPC->Clone(Form("h_non_decay_photons_etaTPC_bin_%02d",i)) );
    vec_non_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_non_decay_photons_etaEMCal->Clone(Form("h_non_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_non_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_non_decay_photons_etaPHOS->Clone(Form("h_non_decay_photons_etaPHOS_bin_%02d",i)) );

    // isolated photon histos TPC
    vec_iso_charged2GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaTPC->Clone(Form("h_iso_charged2GeV_R03_photons_etaTPC_bin_%02d",i)) );
    vec_iso_charged2GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaTPC->Clone(Form("h_iso_charged2GeV_R04_photons_etaTPC_bin_%02d",i)) );
    vec_iso_charged2GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaTPC->Clone(Form("h_iso_charged2GeV_R05_photons_etaTPC_bin_%02d",i)) );
    vec_iso_charged3GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaTPC->Clone(Form("h_iso_charged3GeV_R03_photons_etaTPC_bin_%02d",i)) );
    vec_iso_charged3GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaTPC->Clone(Form("h_iso_charged3GeV_R04_photons_etaTPC_bin_%02d",i)) );
    vec_iso_charged3GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaTPC->Clone(Form("h_iso_charged3GeV_R05_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full2GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaTPC->Clone(Form("h_iso_full2GeV_R03_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full2GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaTPC->Clone(Form("h_iso_full2GeV_R04_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full2GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaTPC->Clone(Form("h_iso_full2GeV_R05_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full3GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaTPC->Clone(Form("h_iso_full3GeV_R03_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full3GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaTPC->Clone(Form("h_iso_full3GeV_R04_photons_etaTPC_bin_%02d",i)) );
    vec_iso_full3GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaTPC->Clone(Form("h_iso_full3GeV_R05_photons_etaTPC_bin_%02d",i)) );

    // isolated photon histos EMCAL
    vec_iso_charged2GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaEMCal->Clone(Form("h_iso_charged2GeV_R03_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_charged2GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaEMCal->Clone(Form("h_iso_charged2GeV_R04_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_charged2GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaEMCal->Clone(Form("h_iso_charged2GeV_R05_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_charged3GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaEMCal->Clone(Form("h_iso_charged3GeV_R03_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_charged3GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaEMCal->Clone(Form("h_iso_charged3GeV_R04_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_charged3GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaEMCal->Clone(Form("h_iso_charged3GeV_R05_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full2GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaEMCal->Clone(Form("h_iso_full2GeV_R03_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full2GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaEMCal->Clone(Form("h_iso_full2GeV_R04_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full2GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaEMCal->Clone(Form("h_iso_full2GeV_R05_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full3GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaEMCal->Clone(Form("h_iso_full3GeV_R03_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full3GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaEMCal->Clone(Form("h_iso_full3GeV_R04_photons_etaEMCal_bin_%02d",i)) );
    vec_iso_full3GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaEMCal->Clone(Form("h_iso_full3GeV_R05_photons_etaEMCal_bin_%02d",i)) );

    // isolated photon histos PHOS
    vec_iso_charged2GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaPHOS->Clone(Form("h_iso_charged2GeV_R03_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_charged2GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaPHOS->Clone(Form("h_iso_charged2GeV_R04_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_charged2GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaPHOS->Clone(Form("h_iso_charged2GeV_R05_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_charged3GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaPHOS->Clone(Form("h_iso_charged3GeV_R03_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_charged3GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaPHOS->Clone(Form("h_iso_charged3GeV_R04_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_charged3GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaPHOS->Clone(Form("h_iso_charged3GeV_R05_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full2GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaPHOS->Clone(Form("h_iso_full2GeV_R03_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full2GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaPHOS->Clone(Form("h_iso_full2GeV_R04_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full2GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaPHOS->Clone(Form("h_iso_full2GeV_R05_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full3GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaPHOS->Clone(Form("h_iso_full3GeV_R03_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full3GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaPHOS->Clone(Form("h_iso_full3GeV_R04_photons_etaPHOS_bin_%02d",i)) );
    vec_iso_full3GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaPHOS->Clone(Form("h_iso_full3GeV_R05_photons_etaPHOS_bin_%02d",i)) );

    // decay photon histos
    vec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(Form("h_decay_photons_etaTPC_bin_%02d",i)) );
    vec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(Form("h_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(Form("h_decay_photons_etaPHOS_bin_%02d",i)) );
    vec_pTHat_bin.push_back( (TH1D*)h_pTHat->Clone(Form("h_pTHat_bin_%02d",i)) );

  }

  //--- begin pTHat bin loop ----------------------------------
  for (int iBin = 0; iBin < pTHatBins; ++iBin) {

    pyHelp.SoftQCD_HardQCD_Switch(iBin, pTHatBin, argv, p, nEvent);
    p.init();

    //--- begin event loop ----------------------------------------------
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      // Generate event.
      if (!p.next()) continue;

      // boost if pPb 5.02 TeV
      if( yesBoost ) p.event.bst(0., 0., 0.435);
      if(iEvent == 0)
	cout << "energy of beam a = " << p.event[1].e() << endl
	     << "energy of beam b = " << p.event[2].e() << endl;

      // reject non-diffractive softQCD events in the hardQCD regime
      if (iBin == 0 && p.info.pTHat() > pTHatBin[1]) continue;

      // reject non-diffractive softQCD events with super large weight, i.e. pthat << ptgamma
      bool is_large_weight = false;
      if(iBin == 0)
	for (int i = 5; i < p.event.size(); i++) {
	  if(p.event[i].isFinal() && p.event[i].id() == 22 && TMath::Abs(p.event[i].eta()) < etaTPC ) {
	    //if(p.event[i].status() < 90 && p.event[i].pT() > 10.)
	    vec_fluctCut.at(0).Fill(p.event[i].pT());
	    for(unsigned int j = 1; j < vec_fluctCut.size(); j++)
	      if(p.event[i].pT() < p.info.pTHat()*(0.9+j*0.1))
		vec_fluctCut.at(j).Fill(p.event[i].pT());
	    
	    if(p.event[i].pT() > p.info.pTHat()*1.7){
		//cout << "softQCD event vetoed with " << "photon pt = " << p.event[i].pT() << "\t and pthat = " << p.info.pTHat() << endl;
	      is_large_weight = true;
	    }
	  }
	}
      if(is_large_weight) continue;

      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaTPC, vec_non_decay_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaEMCal, vec_non_decay_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaPHOS, vec_non_decay_photons_etaPHOS_bin.at(iBin));

      // fill isolated photons: considers only non-decay photons
      // arguments = (p.event, etaAcc, vec_histo, bool onlyCharged?, iso cone radius, iso pt)
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);

      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 3.);


      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);

      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 3.);
      pyHelp.Fill_Non_Decay_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 3.);

      if( !strcmp(argv[2],"decay") ){
        pyHelp.Fill_Decay_Photon_Pt(p.event, etaTPC, vec_decay_photons_etaTPC_bin.at(iBin));
        pyHelp.Fill_Decay_Photon_Pt(p.event, etaEMCal, vec_decay_photons_etaEMCal_bin.at(iBin));
        pyHelp.Fill_Decay_Photon_Pt(p.event, etaPHOS, vec_decay_photons_etaPHOS_bin.at(iBin));
      }

      vec_pTHat_bin.at(iBin)->Fill(p.info.pTHat());

    }// end of event loop

    p.stat();

    double sigma = p.info.sigmaGen()*1e9; // cross section in picobarn
    //    double sigma_per_event = sigma/p.info.weightSum(); // weightSum = number of events in standard Pythia8

    vec_weightSum_bin.at(iBin)->SetBinContent(1,p.info.weightSum());
    h_weightSum->SetBinContent(1,h_weightSum->GetBinContent(1)+p.info.weightSum());
    cout << "- - - weightSum() = " << p.info.weightSum() << endl;

    vec_non_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_non_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_non_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_iso_charged2GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_charged2GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_iso_charged3GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_charged3GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_iso_full2GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_full2GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_iso_full3GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_iso_full3GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    vec_pTHat_bin.at(iBin)->Scale(sigma);

  }// end of pTHat bin loop

  //--- write to root file ---------------------------------------
  TFile file(rootFileName, "RECREATE");

  for( unsigned int i = 0; i < vec_fluctCut.size(); i++){
    vec_fluctCut.at(i).Write();
  }


  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaTPC_bin, h_non_decay_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaEMCal_bin, h_non_decay_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaPHOS_bin, h_non_decay_photons_etaPHOS, file, 2*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaTPC_bin, h_iso_charged2GeV_R03_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaTPC_bin, h_iso_charged2GeV_R04_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaTPC_bin, h_iso_charged2GeV_R05_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaEMCal_bin, h_iso_charged2GeV_R03_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaEMCal_bin, h_iso_charged2GeV_R04_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaEMCal_bin, h_iso_charged2GeV_R05_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaPHOS_bin, h_iso_charged2GeV_R03_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaPHOS_bin, h_iso_charged2GeV_R04_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaPHOS_bin, h_iso_charged2GeV_R05_photons_etaPHOS, file, 2*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaTPC_bin, h_iso_charged3GeV_R03_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaTPC_bin, h_iso_charged3GeV_R04_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaTPC_bin, h_iso_charged3GeV_R05_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaEMCal_bin, h_iso_charged3GeV_R03_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaEMCal_bin, h_iso_charged3GeV_R04_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaEMCal_bin, h_iso_charged3GeV_R05_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaPHOS_bin, h_iso_charged3GeV_R03_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaPHOS_bin, h_iso_charged3GeV_R04_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaPHOS_bin, h_iso_charged3GeV_R05_photons_etaPHOS, file, 2*etaPHOS);


  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaTPC_bin, h_iso_full2GeV_R03_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaTPC_bin, h_iso_full2GeV_R04_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaTPC_bin, h_iso_full2GeV_R05_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaEMCal_bin, h_iso_full2GeV_R03_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaEMCal_bin, h_iso_full2GeV_R04_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaEMCal_bin, h_iso_full2GeV_R05_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaPHOS_bin, h_iso_full2GeV_R03_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaPHOS_bin, h_iso_full2GeV_R04_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaPHOS_bin, h_iso_full2GeV_R05_photons_etaPHOS, file, 2*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaTPC_bin, h_iso_full3GeV_R03_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaTPC_bin, h_iso_full3GeV_R04_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaTPC_bin, h_iso_full3GeV_R05_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaEMCal_bin, h_iso_full3GeV_R03_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaEMCal_bin, h_iso_full3GeV_R04_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaEMCal_bin, h_iso_full3GeV_R05_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaPHOS_bin, h_iso_full3GeV_R03_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaPHOS_bin, h_iso_full3GeV_R04_photons_etaPHOS, file, 2*etaPHOS);
  pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaPHOS_bin, h_iso_full3GeV_R05_photons_etaPHOS, file, 2*etaPHOS);

  if( !strcmp(argv[2],"decay") ){
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaTPC_bin, h_decay_photons_etaTPC, file, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaEMCal_bin, h_decay_photons_etaEMCal, file, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaPHOS_bin, h_decay_photons_etaPHOS, file, 2*etaPHOS);
  }

  pyHelp.Add_Histos_Scale_Write2File( vec_pTHat_bin, h_pTHat, file, 1.);

  for(int iBin=0; iBin < pTHatBins; iBin++){
    vec_weightSum_bin.at(iBin)->Write();
  }
  h_weightSum->Write();

  file.Close();

  return 0;
}
