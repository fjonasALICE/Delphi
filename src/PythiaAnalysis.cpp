#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "PythiaAnalysisHelper.h"
#include "PythiaAnalysis.h"

using std::cout;
using namespace Pythia8;

int main(int, char **);
int main(int argc, char **argv) {

  // not so much output at the beginning
  p.readString("Next:numberCount = 100000");
  p.readString("Next:numberShowLHA = 0");
  p.readString("Next:numberShowInfo = 0");
  p.readString("Next:numberShowProcess = 0");
  p.readString("Next:numberShowEvent = 0");
  pyHelp.Set_Pythia_Randomseed(p);

  //--- read commandline args ----------------------------------------
  if (argc < 5) {
    printf("Need at least first 4 arguments:\n%s [output root file] [\"MB\",\"MBVeto\",\"JJ\",\"GJ\",\"WeakBoson\"] [number of events per pthatbin] [cm energy in GeV] [\"fullEvents\",\"noMPI\",\"noHadro\",\"noMPInoHadro\",\"noShower\"] [renormScaleFac] [factorMultFac] [beta_boost_z] [pdfA] [pdfB]", argv[0]);
    exit(EXIT_FAILURE);
  }

  //--- define output root file --------------------------------------
  // argv[1]: "abc" in output "abc.root"
  // argv[2]: process switch
  // argv[3]: number of events
  // argv[4]: eCM
  // argv[5]: optional argument to switch off MPI, hadronization or entire shower
  // argv[6]: renormMultFac
  // argv[7]: factorMultFac
  // argv[8]: boost in z direction (beta=v/c)
  // argv[9]: external pdf beam A 
  // argv[10]: external pdf beam B

  snprintf( rootFileName, sizeof(rootFileName), "%s.root", argv[1]);  // "abc" -> "abc.root"
  printf("\nThe result will be written into %s\n", rootFileName);

  // p.readString("TimeShower:pTminChgQ = 2.0");  // test shower cut-off

  if (argc >= 9){ // account for boost in asymmetric collision systems
    applyBoost = true;
    boostBetaZ = strtof(argv[8], NULL);
    printf("\nApplying a boost along z direction with beta_z = %f\n", boostBetaZ); // 0.435 for pPb
  }

  if (argc >= 10){ // choice of external PDF (using LHAPDF6)
    string pdfA = argv[9];
    string pdfB = argv[9];
    printf("\nUsing PDF %s for beam A\n", pdfA.c_str());
    p.readString("PDF:pSet = LHAPDF6:" + pdfA);
    if( argc >= 11)
      pdfB = argv[10];
    printf("Using PDF %s for beam B\n", pdfB.c_str());
    p.readString("PDF:pSetB = LHAPDF6:" + pdfB);
  }

  int nEvent = strtol(argv[3], NULL, 10); // number of events
  printf("\nGenerating %d events per pthat bin\n", nEvent);
  
  pyHelp.Pass_Parameters_To_Pythia(p, argc, argv); // which energy, scales, optional master switches

  if (!strcmp(argv[2],"MBVeto")) // prevent double sampling of MB events and e.g. JJ events
    MB_veto = true;

  if( !strcmp(argv[2],"JJ") || !strcmp(argv[2],"GJ") || !strcmp(argv[2],"WeakBoson") ){
    printf("\nUsing %d pTHat bins:\n", pTHatBins);
    for(int i=0; i <= pTHatBins; i++){
      printf("%.0f ", pTHatBin[i]);
    }
    printf("\n");
  }

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------
  //--- Histograms ---------------------------------------------------
  TH1::SetDefaultSumw2(kTRUE); // histograms shall carry sum of weights by default for proper error propagation
  gStyle->SetOptStat(0); // no legend by default

  // TH2D electron_pt vs electron_topMotherID
  TH2D *h2_electron_pt_topMotherID = new TH2D("h2_electron_pt_topMotherID","electron_pt_topMotherID (EMCal acceptance |#eta| < 0.66)",17,0,17,ptBins,ptMin,ptMax);
  h2_electron_pt_topMotherID->SetCanExtend(TH1::kXaxis);
  for(int i = 0; i < 17; i++)
    h2_electron_pt_topMotherID->GetXaxis()->SetBinLabel(i+1, electronMotherName[i]);

  // all pions without secondary correction
  TH1D *h_pi0_etaLarge = new TH1D("h_pi0_etaLarge","pi0_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaTPC   = new TH1D("h_pi0_etaTPC","pi0_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaEMCal = new TH1D("h_pi0_etaEMCal","pi0_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaPHOS  = new TH1D("h_pi0_etaPHOS","pi0_etaPHOS", ptBins, ptMin, ptMax);

  // primary pions (with secondary correction)
  TH1D *h_pi0primary_etaLarge = new TH1D("h_pi0primary_etaLarge","pi0primary_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaTPC   = new TH1D("h_pi0primary_etaTPC","pi0primary_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaEMCal = new TH1D("h_pi0primary_etaEMCal","pi0primary_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaPHOS  = new TH1D("h_pi0primary_etaPHOS","pi0primary_etaPHOS", ptBins, ptMin, ptMax);

  // eta meson
  TH1D *h_eta_etaLarge = new TH1D("h_eta_etaLarge","eta_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaTPC   = new TH1D("h_eta_etaTPC","eta_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaEMCal = new TH1D("h_eta_etaEMCal","eta_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaPHOS  = new TH1D("h_eta_etaPHOS","eta_etaPHOS", ptBins, ptMin, ptMax);

  // eta prime meson
  TH1D *h_etaprime_etaLarge = new TH1D("h_etaprime_etaLarge","etaprime_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaTPC   = new TH1D("h_etaprime_etaTPC","etaprime_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaEMCal = new TH1D("h_etaprime_etaEMCal","etaprime_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaPHOS  = new TH1D("h_etaprime_etaPHOS","etaprime_etaPHOS", ptBins, ptMin, ptMax);

  // omega meson
  TH1D *h_omega_etaLarge = new TH1D("h_omega_etaLarge","omega_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaTPC   = new TH1D("h_omega_etaTPC","omega_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaEMCal = new TH1D("h_omega_etaEMCal","omega_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaPHOS  = new TH1D("h_omega_etaPHOS","omega_etaPHOS", ptBins, ptMin, ptMax);

  // direct photons (consider only direct photons)
  TH1D *h_direct_photons_etaLarge = new TH1D("h_direct_photons_etaLarge","direct photons_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaTPC   = new TH1D("h_direct_photons_etaTPC","direct photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaEMCal = new TH1D("h_direct_photons_etaEMCal","direct photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaPHOS  = new TH1D("h_direct_photons_etaPHOS","direct photons_etaPHOS", ptBins, ptMin, ptMax);

  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_shower_photons_etaLarge = new TH1D("h_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaTPC   = new TH1D("h_shower_photons_etaTPC","shower photons (q -> q #gamma) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaEMCal = new TH1D("h_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaPHOS  = new TH1D("h_shower_photons_etaPHOS","shower photons (q -> q #gamma) in |#eta| < 0.12", ptBins, ptMin, ptMax);
  // photons from ME (aka prompt)
  TH1D *h_222_photons_etaLarge = new TH1D("h_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaTPC   = new TH1D("h_222_photons_etaTPC","photons from ME (aka prompt) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaEMCal = new TH1D("h_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaPHOS  = new TH1D("h_222_photons_etaPHOS","photons from ME (aka prompt) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // isolated photons (considers only direct photons)
  TH1D *h_iso_charged2GeV_R03_photons_etaTPC = new TH1D("h_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R04_photons_etaTPC = new TH1D("h_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R05_photons_etaTPC = new TH1D("h_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaPHOS = new TH1D("h_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R03_photons_etaTPC = new TH1D("h_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R04_photons_etaTPC = new TH1D("h_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R05_photons_etaTPC = new TH1D("h_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaPHOS = new TH1D("h_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R03_photons_etaTPC = new TH1D("h_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaEMCal = new TH1D("h_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaPHOS = new TH1D("h_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R04_photons_etaTPC = new TH1D("h_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaEMCal = new TH1D("h_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaPHOS = new TH1D("h_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R05_photons_etaTPC = new TH1D("h_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaEMCal = new TH1D("h_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaPHOS = new TH1D("h_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R03_photons_etaTPC = new TH1D("h_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaEMCal = new TH1D("h_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaPHOS = new TH1D("h_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R04_photons_etaTPC = new TH1D("h_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaEMCal = new TH1D("h_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaPHOS = new TH1D("h_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R05_photons_etaTPC = new TH1D("h_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaEMCal = new TH1D("h_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaPHOS = new TH1D("h_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  // decay photons
  TH1D *h_decay_photons_etaLarge = new TH1D("h_decay_photons_etaLarge","decay photons_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaTPC = new TH1D("h_decay_photons_etaTPC","decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaEMCal = new TH1D("h_decay_photons_etaEMCal","decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaPHOS = new TH1D("h_decay_photons_etaPHOS","decay photons_etaPHOS", ptBins, ptMin, ptMax);


  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------
  // all pions without secondary correction
  TH1D *h_invXsec_pi0_etaLarge = new TH1D("h_invXsec_pi0_etaLarge","pi0_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaTPC = new TH1D("h_invXsec_pi0_etaTPC","pi0_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaEMCal = new TH1D("h_invXsec_pi0_etaEMCal","pi0_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaPHOS = new TH1D("h_invXsec_pi0_etaPHOS","pi0_etaPHOS", ptBins, ptMin, ptMax);

  // primary pions (with secondary correction)
  TH1D *h_invXsec_pi0primary_etaLarge = new TH1D("h_invXsec_pi0primary_etaLarge","pi0primary_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaTPC = new TH1D("h_invXsec_pi0primary_etaTPC","pi0primary_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaEMCal = new TH1D("h_invXsec_pi0primary_etaEMCal","pi0primary_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaPHOS = new TH1D("h_invXsec_pi0primary_etaPHOS","pi0primary_etaPHOS", ptBins, ptMin, ptMax);

  // eta meson
  TH1D *h_invXsec_eta_etaLarge = new TH1D("h_invXsec_eta_etaLarge","eta_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaTPC = new TH1D("h_invXsec_eta_etaTPC","eta_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaEMCal = new TH1D("h_invXsec_eta_etaEMCal","eta_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaPHOS = new TH1D("h_invXsec_eta_etaPHOS","eta_etaPHOS", ptBins, ptMin, ptMax);

  // eta prime meson
  TH1D *h_invXsec_etaprime_etaLarge = new TH1D("h_invXsec_etaprime_etaLarge","etaprime_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaTPC = new TH1D("h_invXsec_etaprime_etaTPC","etaprime_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaEMCal = new TH1D("h_invXsec_etaprime_etaEMCal","etaprime_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaPHOS = new TH1D("h_invXsec_etaprime_etaPHOS","etaprime_etaPHOS", ptBins, ptMin, ptMax);

  // omega meson
  TH1D *h_invXsec_omega_etaLarge = new TH1D("h_invXsec_omega_etaLarge","omega_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaTPC = new TH1D("h_invXsec_omega_etaTPC","omega_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaEMCal = new TH1D("h_invXsec_omega_etaEMCal","omega_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaPHOS = new TH1D("h_invXsec_omega_etaPHOS","omega_etaPHOS", ptBins, ptMin, ptMax);

  // direct photons (consider only direct photons)
  TH1D *h_invXsec_direct_photons_etaLarge = new TH1D("h_invXsec_direct_photons_etaLarge","direct photons_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaTPC = new TH1D("h_invXsec_direct_photons_etaTPC","direct photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaEMCal = new TH1D("h_invXsec_direct_photons_etaEMCal","direct photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaPHOS = new TH1D("h_invXsec_direct_photons_etaPHOS","direct photons_etaPHOS", ptBins, ptMin, ptMax);

  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_invXsec_shower_photons_etaLarge = new TH1D("h_invXsec_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaTPC   = new TH1D("h_invXsec_shower_photons_etaTPC","shower photons (q -> q #gamma) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaEMCal = new TH1D("h_invXsec_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaPHOS  = new TH1D("h_invXsec_shower_photons_etaPHOS","shower photons (q -> q #gamma) in |#eta| < 0.12", ptBins, ptMin, ptMax);
  // photons from ME (aka prompt)
  TH1D *h_invXsec_222_photons_etaLarge = new TH1D("h_invXsec_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaTPC   = new TH1D("h_invXsec_222_photons_etaTPC","photons from ME (aka prompt) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaEMCal = new TH1D("h_invXsec_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaPHOS  = new TH1D("h_invXsec_222_photons_etaPHOS","photons from ME (aka prompt) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // isolated photons (considers only direct photons)
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaTPC = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaPHOS = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaTPC = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaPHOS = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaTPC = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaPHOS = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaTPC = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaPHOS = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaTPC = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaPHOS = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaTPC = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaPHOS = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaTPC = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaPHOS = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaTPC = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaPHOS = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaTPC = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaPHOS = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaTPC = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaPHOS = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaTPC = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaPHOS = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaTPC = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaPHOS = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);
  
  // decay photons
  TH1D *h_invXsec_decay_photons_etaLarge = new TH1D("h_invXsec_decay_photons_etaLarge","decay photons_etaLarge", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaTPC = new TH1D("h_invXsec_decay_photons_etaTPC","decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaEMCal = new TH1D("h_invXsec_decay_photons_etaEMCal","decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaPHOS = new TH1D("h_invXsec_decay_photons_etaPHOS","decay photons_etaPHOS", ptBins, ptMin, ptMax);



  //----------------------------------------------------------------------------------------------------
  // check underlying born kt to see, e.g., if HardQCD cross section does not blow up
  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", 1000, ptMin, ptMax);
  // store the weight sum for proper normalization afterwards
  TH1D *h_weightSum = new TH1D("h_weightSum","sum of weights", 1, 0., 1.);

  // organise pTHat wise histograms in vectors
  vector <TH2D*> vec_electron_pt_topMotherID_bin;

  vector <TH1D*> vec_pi0_etaLarge_bin;
  vector <TH1D*> vec_pi0_etaTPC_bin;
  vector <TH1D*> vec_pi0_etaEMCal_bin;
  vector <TH1D*> vec_pi0_etaPHOS_bin;

  vector <TH1D*> vec_pi0primary_etaLarge_bin;
  vector <TH1D*> vec_pi0primary_etaTPC_bin;
  vector <TH1D*> vec_pi0primary_etaEMCal_bin;
  vector <TH1D*> vec_pi0primary_etaPHOS_bin;

  vector <TH1D*> vec_eta_etaLarge_bin;
  vector <TH1D*> vec_eta_etaTPC_bin;
  vector <TH1D*> vec_eta_etaEMCal_bin;
  vector <TH1D*> vec_eta_etaPHOS_bin;

  vector <TH1D*> vec_etaprime_etaLarge_bin;
  vector <TH1D*> vec_etaprime_etaTPC_bin;
  vector <TH1D*> vec_etaprime_etaEMCal_bin;
  vector <TH1D*> vec_etaprime_etaPHOS_bin;

  vector <TH1D*> vec_omega_etaLarge_bin;
  vector <TH1D*> vec_omega_etaTPC_bin;
  vector <TH1D*> vec_omega_etaEMCal_bin;
  vector <TH1D*> vec_omega_etaPHOS_bin;

  vector <TH1D*> vec_direct_photons_etaLarge_bin;
  vector <TH1D*> vec_direct_photons_etaTPC_bin;
  vector <TH1D*> vec_direct_photons_etaEMCal_bin;
  vector <TH1D*> vec_direct_photons_etaPHOS_bin;

  vector <TH1D*> vec_shower_photons_etaLarge_bin;
  vector <TH1D*> vec_shower_photons_etaTPC_bin;
  vector <TH1D*> vec_shower_photons_etaEMCal_bin;
  vector <TH1D*> vec_shower_photons_etaPHOS_bin;

  vector <TH1D*> vec_222_photons_etaLarge_bin;
  vector <TH1D*> vec_222_photons_etaTPC_bin;
  vector <TH1D*> vec_222_photons_etaEMCal_bin;
  vector <TH1D*> vec_222_photons_etaPHOS_bin;

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

  vector <TH1D*> vec_decay_photons_etaLarge_bin;
  vector <TH1D*> vec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_decay_photons_etaPHOS_bin;

  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------
  vector <TH1D*> vec_invXsec_pi0_etaLarge_bin;
  vector <TH1D*> vec_invXsec_pi0_etaTPC_bin;
  vector <TH1D*> vec_invXsec_pi0_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_pi0_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_pi0primary_etaLarge_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaTPC_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_eta_etaLarge_bin;
  vector <TH1D*> vec_invXsec_eta_etaTPC_bin;
  vector <TH1D*> vec_invXsec_eta_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_eta_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_etaprime_etaLarge_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaTPC_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_omega_etaLarge_bin;
  vector <TH1D*> vec_invXsec_omega_etaTPC_bin;
  vector <TH1D*> vec_invXsec_omega_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_omega_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_direct_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_shower_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_222_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin;

  vector <TH1D*> vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin;

  vector <TH1D*> vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin;
  vector <TH1D*> vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_decay_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaPHOS_bin;

  //----------------------------------------------------------------------------------------------------

  vector <TH1D*> vec_pTHat_bin;
  vector <TH1D*> vec_weightSum_bin;

  for(int i = 0; i < pTHatBins; i++){

    vec_electron_pt_topMotherID_bin.push_back( (TH2D*)h2_electron_pt_topMotherID->Clone(Form( "h2_electron_pt_topMotherID_bin_%02d", i)) );

    vec_weightSum_bin.push_back( (TH1D*)h_weightSum->Clone(Form( "h_weightSum_bin_%02d", i )) );

    // pi0 histos
    vec_pi0_etaLarge_bin.push_back( (TH1D*)h_pi0_etaLarge->Clone(Form("h_pi0_etaLarge_bin_%02d",i)) );
    vec_pi0_etaTPC_bin.push_back( (TH1D*)h_pi0_etaTPC->Clone(Form("h_pi0_etaTPC_bin_%02d",i)) );
    vec_pi0_etaEMCal_bin.push_back( (TH1D*)h_pi0_etaEMCal->Clone(Form("h_pi0_etaEMCal_bin_%02d",i)) );
    vec_pi0_etaPHOS_bin.push_back( (TH1D*)h_pi0_etaPHOS->Clone(Form("h_pi0_etaPHOS_bin_%02d",i)) );

    // primary pi0 histos
    vec_pi0primary_etaLarge_bin.push_back( (TH1D*)h_pi0primary_etaLarge->Clone(Form("h_pi0primary_etaLarge_bin_%02d",i)) );
    vec_pi0primary_etaTPC_bin.push_back( (TH1D*)h_pi0primary_etaTPC->Clone(Form("h_pi0primary_etaTPC_bin_%02d",i)) );
    vec_pi0primary_etaEMCal_bin.push_back( (TH1D*)h_pi0primary_etaEMCal->Clone(Form("h_pi0primary_etaEMCal_bin_%02d",i)) );
    vec_pi0primary_etaPHOS_bin.push_back( (TH1D*)h_pi0primary_etaPHOS->Clone(Form("h_pi0primary_etaPHOS_bin_%02d",i)) );

    // eta histos
    vec_eta_etaLarge_bin.push_back( (TH1D*)h_eta_etaLarge->Clone(Form("h_eta_etaLarge_bin_%02d",i)) );
    vec_eta_etaTPC_bin.push_back( (TH1D*)h_eta_etaTPC->Clone(Form("h_eta_etaTPC_bin_%02d",i)) );
    vec_eta_etaEMCal_bin.push_back( (TH1D*)h_eta_etaEMCal->Clone(Form("h_eta_etaEMCal_bin_%02d",i)) );
    vec_eta_etaPHOS_bin.push_back( (TH1D*)h_eta_etaPHOS->Clone(Form("h_eta_etaPHOS_bin_%02d",i)) );

    // eta prime histos
    vec_etaprime_etaLarge_bin.push_back( (TH1D*)h_etaprime_etaLarge->Clone(Form("h_etaprime_etaLarge_bin_%02d",i)) );
    vec_etaprime_etaTPC_bin.push_back( (TH1D*)h_etaprime_etaTPC->Clone(Form("h_etaprime_etaTPC_bin_%02d",i)) );
    vec_etaprime_etaEMCal_bin.push_back( (TH1D*)h_etaprime_etaEMCal->Clone(Form("h_etaprime_etaEMCal_bin_%02d",i)) );
    vec_etaprime_etaPHOS_bin.push_back( (TH1D*)h_etaprime_etaPHOS->Clone(Form("h_etaprime_etaPHOS_bin_%02d",i)) );

    // omega histos
    vec_omega_etaLarge_bin.push_back( (TH1D*)h_omega_etaLarge->Clone(Form("h_omega_etaLarge_bin_%02d",i)) );
    vec_omega_etaTPC_bin.push_back( (TH1D*)h_omega_etaTPC->Clone(Form("h_omega_etaTPC_bin_%02d",i)) );
    vec_omega_etaEMCal_bin.push_back( (TH1D*)h_omega_etaEMCal->Clone(Form("h_omega_etaEMCal_bin_%02d",i)) );
    vec_omega_etaPHOS_bin.push_back( (TH1D*)h_omega_etaPHOS->Clone(Form("h_omega_etaPHOS_bin_%02d",i)) );

    // direct photon histos
    vec_direct_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_direct_photons_etaLarge_bin_%02d",i)) );
    vec_direct_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_direct_photons_etaTPC_bin_%02d",i)) );
    vec_direct_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_direct_photons_etaEMCal_bin_%02d",i)) );
    vec_direct_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_direct_photons_etaPHOS_bin_%02d",i)) );

    // discriminating shower and prompt photons
    vec_shower_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_shower_photons_etaLarge_bin_%02d",i)) );
    vec_shower_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_shower_photons_etaTPC_bin_%02d",i)) );
    vec_shower_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_shower_photons_etaEMCal_bin_%02d",i)) );
    vec_shower_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_shower_photons_etaPHOS_bin_%02d",i)) );

    vec_222_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_222_photons_etaLarge_bin_%02d",i)) );
    vec_222_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_222_photons_etaTPC_bin_%02d",i)) );
    vec_222_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_222_photons_etaEMCal_bin_%02d",i)) );
    vec_222_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_222_photons_etaPHOS_bin_%02d",i)) );

    if(applyPhotonIso){
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
    }

    // decay photon histos
    vec_decay_photons_etaLarge_bin.push_back( (TH1D*)h_decay_photons_etaLarge->Clone(Form("h_decay_photons_etaLarge_bin_%02d",i)) );
    vec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(Form("h_decay_photons_etaTPC_bin_%02d",i)) );
    vec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(Form("h_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(Form("h_decay_photons_etaPHOS_bin_%02d",i)) );



    //----------------------------------------------------------------------------------------------------
    // do the same jazz for invariant cross section histos ------------------------------------------
    //------------------------------------------------------------------------------------------
    // pi0 histos
    vec_invXsec_pi0_etaLarge_bin.push_back( (TH1D*)h_invXsec_pi0_etaLarge->Clone(Form("h_invXsec_pi0_etaLarge_bin_%02d",i)) );
    vec_invXsec_pi0_etaTPC_bin.push_back( (TH1D*)h_invXsec_pi0_etaTPC->Clone(Form("h_invXsec_pi0_etaTPC_bin_%02d",i)) );
    vec_invXsec_pi0_etaEMCal_bin.push_back( (TH1D*)h_invXsec_pi0_etaEMCal->Clone(Form("h_invXsec_pi0_etaEMCal_bin_%02d",i)) );
    vec_invXsec_pi0_etaPHOS_bin.push_back( (TH1D*)h_invXsec_pi0_etaPHOS->Clone(Form("h_invXsec_pi0_etaPHOS_bin_%02d",i)) );

    // primary pi0 histos
    vec_invXsec_pi0primary_etaLarge_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaLarge->Clone(Form("h_invXsec_pi0primary_etaLarge_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaTPC_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaTPC->Clone(Form("h_invXsec_pi0primary_etaTPC_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaEMCal_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaEMCal->Clone(Form("h_invXsec_pi0primary_etaEMCal_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaPHOS_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaPHOS->Clone(Form("h_invXsec_pi0primary_etaPHOS_bin_%02d",i)) );

    // eta histos
    vec_invXsec_eta_etaLarge_bin.push_back( (TH1D*)h_invXsec_eta_etaLarge->Clone(Form("h_invXsec_eta_etaLarge_bin_%02d",i)) );
    vec_invXsec_eta_etaTPC_bin.push_back( (TH1D*)h_invXsec_eta_etaTPC->Clone(Form("h_invXsec_eta_etaTPC_bin_%02d",i)) );
    vec_invXsec_eta_etaEMCal_bin.push_back( (TH1D*)h_invXsec_eta_etaEMCal->Clone(Form("h_invXsec_eta_etaEMCal_bin_%02d",i)) );
    vec_invXsec_eta_etaPHOS_bin.push_back( (TH1D*)h_invXsec_eta_etaPHOS->Clone(Form("h_invXsec_eta_etaPHOS_bin_%02d",i)) );

    // eta prime histos
    vec_invXsec_etaprime_etaLarge_bin.push_back( (TH1D*)h_invXsec_etaprime_etaLarge->Clone(Form("h_invXsec_etaprime_etaLarge_bin_%02d",i)) );
    vec_invXsec_etaprime_etaTPC_bin.push_back( (TH1D*)h_invXsec_etaprime_etaTPC->Clone(Form("h_invXsec_etaprime_etaTPC_bin_%02d",i)) );
    vec_invXsec_etaprime_etaEMCal_bin.push_back( (TH1D*)h_invXsec_etaprime_etaEMCal->Clone(Form("h_invXsec_etaprime_etaEMCal_bin_%02d",i)) );
    vec_invXsec_etaprime_etaPHOS_bin.push_back( (TH1D*)h_invXsec_etaprime_etaPHOS->Clone(Form("h_invXsec_etaprime_etaPHOS_bin_%02d",i)) );

    // omega histos
    vec_invXsec_omega_etaLarge_bin.push_back( (TH1D*)h_invXsec_omega_etaLarge->Clone(Form("h_invXsec_omega_etaLarge_bin_%02d",i)) );
    vec_invXsec_omega_etaTPC_bin.push_back( (TH1D*)h_invXsec_omega_etaTPC->Clone(Form("h_invXsec_omega_etaTPC_bin_%02d",i)) );
    vec_invXsec_omega_etaEMCal_bin.push_back( (TH1D*)h_invXsec_omega_etaEMCal->Clone(Form("h_invXsec_omega_etaEMCal_bin_%02d",i)) );
    vec_invXsec_omega_etaPHOS_bin.push_back( (TH1D*)h_invXsec_omega_etaPHOS->Clone(Form("h_invXsec_omega_etaPHOS_bin_%02d",i)) );

    // direct photon histos
    vec_invXsec_direct_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_direct_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_direct_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_direct_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_direct_photons_etaPHOS_bin_%02d",i)) );

    // discriminating shower and prompt photons
    vec_invXsec_shower_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_shower_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_shower_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_shower_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_shower_photons_etaPHOS_bin_%02d",i)) );

    vec_invXsec_222_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_222_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_222_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_222_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_222_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_222_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_222_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_222_photons_etaPHOS_bin_%02d",i)) );

    if(applyPhotonIso){
      // isolated photon histos TPC
      vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaTPC->Clone(Form("h_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaTPC->Clone(Form("h_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaTPC->Clone(Form("h_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaTPC->Clone(Form("h_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaTPC->Clone(Form("h_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaTPC->Clone(Form("h_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaTPC->Clone(Form("h_invXsec_iso_full2GeV_R03_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaTPC->Clone(Form("h_invXsec_iso_full2GeV_R04_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaTPC->Clone(Form("h_invXsec_iso_full2GeV_R05_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaTPC->Clone(Form("h_invXsec_iso_full3GeV_R03_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaTPC->Clone(Form("h_invXsec_iso_full3GeV_R04_photons_etaTPC_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaTPC->Clone(Form("h_invXsec_iso_full3GeV_R05_photons_etaTPC_bin_%02d",i)) );

      // isolated photon histos EMCAL
      vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaEMCal->Clone(Form("h_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaEMCal->Clone(Form("h_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaEMCal->Clone(Form("h_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaEMCal->Clone(Form("h_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaEMCal->Clone(Form("h_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaEMCal->Clone(Form("h_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaEMCal->Clone(Form("h_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin_%02d",i)) );

      // isolated photon histos PHOS
      vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R03_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R04_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged2GeV_R05_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R03_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R04_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_charged3GeV_R05_photons_etaPHOS->Clone(Form("h_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R03_photons_etaPHOS->Clone(Form("h_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R04_photons_etaPHOS->Clone(Form("h_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full2GeV_R05_photons_etaPHOS->Clone(Form("h_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R03_photons_etaPHOS->Clone(Form("h_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R04_photons_etaPHOS->Clone(Form("h_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin_%02d",i)) );
      vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin.push_back( (TH1D*)h_iso_full3GeV_R05_photons_etaPHOS->Clone(Form("h_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin_%02d",i)) );
    }

    // decay photon histos
    vec_invXsec_decay_photons_etaLarge_bin.push_back( (TH1D*)h_decay_photons_etaLarge->Clone(Form("h_invXsec_decay_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(Form("h_invXsec_decay_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(Form("h_invXsec_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(Form("h_invXsec_decay_photons_etaPHOS_bin_%02d",i)) );





    //----------------------------------------------------------------------------------------------------
    vec_pTHat_bin.push_back( (TH1D*)h_pTHat->Clone(Form("h_pTHat_bin_%02d",i)) );

  }

  //--- begin pTHat bin loop ----------------------------------
  for (int iBin = pTHatStartBin; iBin < pTHatBins; ++iBin) {

    pyHelp.ProcessSwitch(iBin, pTHatBin, argv, p);
   
    p.init();

    //--- begin event loop ----------------------------------------------
    for (int iEvent = 1; iEvent <= nEvent; ++iEvent) {
      // Generate event.
      if (!p.next()) continue;

      // boost if pPb
      if( applyBoost ) p.event.bst(0., 0., boostBetaZ);
      if(iEvent == 1)
        cout << "energy of beam a = " << p.event[1].e() << endl
             << "energy of beam b = " << p.event[2].e() << endl;

      if ( !strcmp(argv[2],"MBVeto") && MB_veto ) {	//---------------------------------------------------------
	// reject softQCD events in the hardQCD regime
	if (p.info.pTHat() > pTHatBin[iBin]) continue;
	
	// #### omitted (at least for the moment) because it may change cross section ###########
	/*	// reject softQCD events with super large weight, i.e. pthat << ptparticle
	bool is_large_weight = false;
	for (int i = 5; i < p.event.size(); i++) {
	  if(p.event[i].isFinal() && p.info.pTHat() > 5. ) {
	    if(p.event[i].pT() > p.info.pTHat()*3.0){
	      cout << "softQCD event vetoed with " << "final particle pt = " << p.event[i].pT() << "\t and pthat = " << p.info.pTHat() << endl;
	      is_large_weight = true;
	    }
	  }
	}      
	if(is_large_weight) continue;
	*/
	
      }	//------------------------------------------------------------------------------------------



      pyHelp.Fill_TH2_Electron_TopMotherID(p.event, etaEMCal, vec_electron_pt_topMotherID_bin.at(iBin));

      pyHelp.Fill_Pi0_Pt(p.event, etaLarge, vec_pi0_etaLarge_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaTPC, vec_pi0_etaTPC_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaEMCal, vec_pi0_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaPHOS, vec_pi0_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Pi0Primary_Pt(p.event, etaLarge, vec_pi0primary_etaLarge_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaTPC, vec_pi0primary_etaTPC_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaEMCal, vec_pi0primary_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaPHOS, vec_pi0primary_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Eta_Pt(p.event, etaLarge, vec_eta_etaLarge_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaTPC, vec_eta_etaTPC_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaEMCal, vec_eta_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaPHOS, vec_eta_etaPHOS_bin.at(iBin));

      pyHelp.Fill_EtaPrime_Pt(p.event, etaLarge, vec_etaprime_etaLarge_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaTPC, vec_etaprime_etaTPC_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaEMCal, vec_etaprime_etaEMCal_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaPHOS, vec_etaprime_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Omega_Pt(p.event, etaLarge, vec_omega_etaLarge_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaTPC, vec_omega_etaTPC_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaEMCal, vec_omega_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaPHOS, vec_omega_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Direct_Photon_Pt(p.event, etaLarge, vec_direct_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaTPC, vec_direct_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaEMCal, vec_direct_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaPHOS, vec_direct_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Shower_Photon_Pt(p.event, etaLarge, vec_shower_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaTPC, vec_shower_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaEMCal, vec_shower_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaPHOS, vec_shower_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_222_Photon_Pt(p.event, etaLarge, vec_222_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaTPC, vec_222_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaEMCal, vec_222_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaPHOS, vec_222_photons_etaPHOS_bin.at(iBin));

      if(applyPhotonIso){
	// fill isolated photons: considers only direct photons
	// arguments = (p.event, etaAcc, vec_histo, bool onlyCharged?, iso cone radius, iso pt)
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged2GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged2GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged2GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);

	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_charged3GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_charged3GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_charged3GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 3.);


	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full2GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full2GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full2GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);

	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_iso_full3GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_iso_full3GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_iso_full3GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 3.);
      }

      pyHelp.Fill_Decay_Photon_Pt(p.event, etaLarge, vec_decay_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaTPC, vec_decay_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaEMCal, vec_decay_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaPHOS, vec_decay_photons_etaPHOS_bin.at(iBin));

      //----------------------------------------------------------------------------------------------------
      // do the same jazz for invariant cross section histos ------------------------------------------
      //------------------------------------------------------------------------------------------

      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaLarge, vec_invXsec_pi0_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaTPC, vec_invXsec_pi0_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaEMCal, vec_invXsec_pi0_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaPHOS, vec_invXsec_pi0_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaLarge, vec_invXsec_pi0primary_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaTPC, vec_invXsec_pi0primary_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaEMCal, vec_invXsec_pi0primary_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaPHOS, vec_invXsec_pi0primary_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaLarge, vec_invXsec_eta_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaTPC, vec_invXsec_eta_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaEMCal, vec_invXsec_eta_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaPHOS, vec_invXsec_eta_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaLarge, vec_invXsec_etaprime_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaTPC, vec_invXsec_etaprime_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaEMCal, vec_invXsec_etaprime_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaPHOS, vec_invXsec_etaprime_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaLarge, vec_invXsec_omega_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaTPC, vec_invXsec_omega_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaEMCal, vec_invXsec_omega_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaPHOS, vec_invXsec_omega_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaLarge, vec_invXsec_direct_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaTPC, vec_invXsec_direct_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaEMCal, vec_invXsec_direct_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaPHOS, vec_invXsec_direct_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaLarge, vec_invXsec_shower_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaTPC, vec_invXsec_shower_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaEMCal, vec_invXsec_shower_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaPHOS, vec_invXsec_shower_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaLarge, vec_invXsec_222_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaTPC, vec_invXsec_222_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaEMCal, vec_invXsec_222_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaPHOS, vec_invXsec_222_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Decay_Photon_Pt(p.event, etaLarge, vec_invXsec_decay_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Decay_Photon_Pt(p.event, etaTPC, vec_invXsec_decay_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Decay_Photon_Pt(p.event, etaEMCal, vec_invXsec_decay_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Decay_Photon_Pt(p.event, etaPHOS, vec_invXsec_decay_photons_etaPHOS_bin.at(iBin));

      if(applyPhotonIso){
	// fill isolated photons: considers only direct photons
	// arguments = (p.event, etaAcc, vec_invXsec_histo, bool onlyCharged?, iso cone radius, iso pt)
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);

	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin.at(iBin), true, 0.5, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin.at(iBin), true, 0.5, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin.at(iBin), true, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin.at(iBin), true, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin.at(iBin), true, 0.5, 3.);


	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);

	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaTPC, vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin.at(iBin), false, 0.5, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaEMCal, vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin.at(iBin), false, 0.5, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin.at(iBin), false, 0.3, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin.at(iBin), false, 0.4, 3.);
	pyHelp.Fill_invXsec_Direct_Iso_Photon_Pt(p.event, etaPHOS, vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin.at(iBin), false, 0.5, 3.);
      }
      

      //----------------------------------------------------------------------------------------------------
      vec_pTHat_bin.at(iBin)->Fill(p.info.pTHat());

    }// end of event loop

    p.stat();

    double sigma = p.info.sigmaGen()*1e9; // cross section in picobarn
    //    double sigma_per_event = sigma/p.info.weightSum(); // weightSum = number of events in standard Pythia8

    vec_weightSum_bin.at(iBin)->SetBinContent(1,p.info.weightSum());
    h_weightSum->SetBinContent(1,h_weightSum->GetBinContent(1)+p.info.weightSum());
    cout << "- - - weightSum() = " << p.info.weightSum() << endl;

    vec_electron_pt_topMotherID_bin.at(iBin)->Scale(sigma);

    vec_pi0_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_pi0primary_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_eta_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_eta_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_eta_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_eta_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_etaprime_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_omega_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_omega_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_omega_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_omega_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_direct_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_shower_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_222_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaPHOS_bin.at(iBin)->Scale(sigma);


    if(applyPhotonIso){
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
    }

    vec_decay_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    //----------------------------------------------------------------------------------------------------
    // do the same jazz for invariant cross section histos ------------------------------------------
    //------------------------------------------------------------------------------------------

    vec_invXsec_pi0_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_pi0primary_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_eta_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_etaprime_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_omega_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_direct_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_shower_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_222_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    if(applyPhotonIso){
      vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
    }

    vec_invXsec_decay_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);


    //----------------------------------------------------------------------------------------------------
    vec_pTHat_bin.at(iBin)->Scale(sigma);

    //---consider only first pthat bin for MB, because pthatbins do not apply for MB---
    if ( (!strcmp(argv[2],"MB") || !strcmp(argv[2],"MBVeto"))
	 && iBin >= 0 )
      break;
  }// end of pTHat bin loop

  //--- write to root file ---------------------------------------
  TFile file(rootFileName, "RECREATE");
  
  //----------------------------------------------------------------------------------------------------
  for(int iBin=0; iBin < pTHatBins; iBin++){
    vec_weightSum_bin.at(iBin)->Write();
  }
  h_weightSum->Write();

  TDirectory *dir_electron = file.mkdir("electron");
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_pt_topMotherID_bin, h2_electron_pt_topMotherID, file, dir_electron, 2*etaEMCal);

  TDirectory *dir_pTHat = file.mkdir("pTHat");
  pyHelp.Add_Histos_Scale_Write2File( vec_pTHat_bin, h_pTHat, file, dir_pTHat, 1.);

  TDirectory *dir_pi0 = file.mkdir("pi0");
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaLarge_bin, h_pi0_etaLarge, file, dir_pi0, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaTPC_bin, h_pi0_etaTPC, file, dir_pi0, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaEMCal_bin, h_pi0_etaEMCal, file, dir_pi0, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaPHOS_bin, h_pi0_etaPHOS, file, dir_pi0, 2*etaPHOS);

  TDirectory *dir_pi0primary = file.mkdir("pi0primary");
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaLarge_bin, h_pi0primary_etaLarge, file, dir_pi0primary, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaTPC_bin, h_pi0primary_etaTPC, file, dir_pi0primary, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaEMCal_bin, h_pi0primary_etaEMCal, file, dir_pi0primary, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaPHOS_bin, h_pi0primary_etaPHOS, file, dir_pi0primary, 2*etaPHOS);

  TDirectory *dir_eta = file.mkdir("eta");
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaLarge_bin, h_eta_etaLarge, file, dir_eta, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaTPC_bin, h_eta_etaTPC, file, dir_eta, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaEMCal_bin, h_eta_etaEMCal, file, dir_eta, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaPHOS_bin, h_eta_etaPHOS, file, dir_eta, 2*etaPHOS);

  TDirectory *dir_etaprime = file.mkdir("etaprime");
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaLarge_bin, h_etaprime_etaLarge, file, dir_etaprime, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaTPC_bin, h_etaprime_etaTPC, file, dir_etaprime, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaEMCal_bin, h_etaprime_etaEMCal, file, dir_etaprime, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaPHOS_bin, h_etaprime_etaPHOS, file, dir_etaprime, 2*etaPHOS);

  TDirectory *dir_omega = file.mkdir("omega");
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaLarge_bin, h_omega_etaLarge, file, dir_omega, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaTPC_bin, h_omega_etaTPC, file, dir_omega, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaEMCal_bin, h_omega_etaEMCal, file, dir_omega, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaPHOS_bin, h_omega_etaPHOS, file, dir_omega, 2*etaPHOS);

  TDirectory *dir_gamma = file.mkdir("gamma");
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaLarge_bin, h_direct_photons_etaLarge, file, dir_gamma, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaTPC_bin, h_direct_photons_etaTPC, file, dir_gamma, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaEMCal_bin, h_direct_photons_etaEMCal, file, dir_gamma, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaPHOS_bin, h_direct_photons_etaPHOS, file, dir_gamma, 2*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaLarge_bin, h_shower_photons_etaLarge, file, dir_gamma, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaTPC_bin, h_shower_photons_etaTPC, file, dir_gamma, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaEMCal_bin, h_shower_photons_etaEMCal, file, dir_gamma, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaPHOS_bin, h_shower_photons_etaPHOS, file, dir_gamma, 2*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaLarge_bin, h_222_photons_etaLarge, file, dir_gamma, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaTPC_bin, h_222_photons_etaTPC, file, dir_gamma, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaEMCal_bin, h_222_photons_etaEMCal, file, dir_gamma, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaPHOS_bin, h_222_photons_etaPHOS, file, dir_gamma, 2*etaPHOS);

  if(applyPhotonIso){
    TDirectory *dir_isoGamma = file.mkdir("isoGamma");
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaTPC_bin, h_iso_charged2GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaTPC_bin, h_iso_charged2GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaTPC_bin, h_iso_charged2GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaEMCal_bin, h_iso_charged2GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaEMCal_bin, h_iso_charged2GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaEMCal_bin, h_iso_charged2GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaPHOS_bin, h_iso_charged2GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaPHOS_bin, h_iso_charged2GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaPHOS_bin, h_iso_charged2GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);

    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaTPC_bin, h_iso_charged3GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaTPC_bin, h_iso_charged3GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaTPC_bin, h_iso_charged3GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaEMCal_bin, h_iso_charged3GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaEMCal_bin, h_iso_charged3GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaEMCal_bin, h_iso_charged3GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaPHOS_bin, h_iso_charged3GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaPHOS_bin, h_iso_charged3GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaPHOS_bin, h_iso_charged3GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);


    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaTPC_bin, h_iso_full2GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaTPC_bin, h_iso_full2GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaTPC_bin, h_iso_full2GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaEMCal_bin, h_iso_full2GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaEMCal_bin, h_iso_full2GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaEMCal_bin, h_iso_full2GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaPHOS_bin, h_iso_full2GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaPHOS_bin, h_iso_full2GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaPHOS_bin, h_iso_full2GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);

    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaTPC_bin, h_iso_full3GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaTPC_bin, h_iso_full3GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaTPC_bin, h_iso_full3GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaEMCal_bin, h_iso_full3GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaEMCal_bin, h_iso_full3GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaEMCal_bin, h_iso_full3GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaPHOS_bin, h_iso_full3GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaPHOS_bin, h_iso_full3GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaPHOS_bin, h_iso_full3GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS);
  }

  TDirectory *dir_decayGamma = file.mkdir("decayGamma");
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaLarge_bin, h_decay_photons_etaLarge, file, dir_decayGamma, 2*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaTPC_bin, h_decay_photons_etaTPC, file, dir_decayGamma, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaEMCal_bin, h_decay_photons_etaEMCal, file, dir_decayGamma, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaPHOS_bin, h_decay_photons_etaPHOS, file, dir_decayGamma, 2*etaPHOS);
  

  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------

  TDirectory *dir_pi0_invXsec = file.mkdir("pi0_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaLarge_bin, h_invXsec_pi0_etaLarge, file, dir_pi0_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaTPC_bin, h_invXsec_pi0_etaTPC, file, dir_pi0_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaEMCal_bin, h_invXsec_pi0_etaEMCal, file, dir_pi0_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaPHOS_bin, h_invXsec_pi0_etaPHOS, file, dir_pi0_invXsec, 2.*etaPHOS);

  TDirectory *dir_pi0primary_invXsec = file.mkdir("pi0primary_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaLarge_bin, h_invXsec_pi0primary_etaLarge, file, dir_pi0primary_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaTPC_bin, h_invXsec_pi0primary_etaTPC, file, dir_pi0primary_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaEMCal_bin, h_invXsec_pi0primary_etaEMCal, file, dir_pi0primary_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaPHOS_bin, h_invXsec_pi0primary_etaPHOS, file, dir_pi0primary_invXsec, 2.*etaPHOS);

  TDirectory *dir_eta_invXsec = file.mkdir("eta_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaLarge_bin, h_invXsec_eta_etaLarge, file, dir_eta_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaTPC_bin, h_invXsec_eta_etaTPC, file, dir_eta_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaEMCal_bin, h_invXsec_eta_etaEMCal, file, dir_eta_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaPHOS_bin, h_invXsec_eta_etaPHOS, file, dir_eta_invXsec, 2.*etaPHOS);

  TDirectory *dir_etaprime_invXsec = file.mkdir("etaprime_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaLarge_bin, h_invXsec_etaprime_etaLarge, file, dir_etaprime_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaTPC_bin, h_invXsec_etaprime_etaTPC, file, dir_etaprime_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaEMCal_bin, h_invXsec_etaprime_etaEMCal, file, dir_etaprime_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaPHOS_bin, h_invXsec_etaprime_etaPHOS, file, dir_etaprime_invXsec, 2.*etaPHOS);

  TDirectory *dir_omega_invXsec = file.mkdir("omega_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaLarge_bin, h_invXsec_omega_etaLarge, file, dir_omega_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaTPC_bin, h_invXsec_omega_etaTPC, file, dir_omega_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaEMCal_bin, h_invXsec_omega_etaEMCal, file, dir_omega_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaPHOS_bin, h_invXsec_omega_etaPHOS, file, dir_omega_invXsec, 2.*etaPHOS);

  TDirectory *dir_gamma_invXsec = file.mkdir("gamma_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaLarge_bin, h_invXsec_direct_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaTPC_bin, h_invXsec_direct_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaEMCal_bin, h_invXsec_direct_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaPHOS_bin, h_invXsec_direct_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaLarge_bin, h_invXsec_shower_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaTPC_bin, h_invXsec_shower_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaEMCal_bin, h_invXsec_shower_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaPHOS_bin, h_invXsec_shower_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS);

  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaLarge_bin, h_invXsec_222_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaTPC_bin, h_invXsec_222_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaEMCal_bin, h_invXsec_222_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaPHOS_bin, h_invXsec_222_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS);

  if(applyPhotonIso){
    TDirectory *dir_isoGamma_invXsec = file.mkdir("isoGamma_invXsec");
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);

    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);


    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);

    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS);
  }

  TDirectory *dir_decayGamma_invXsec = file.mkdir("decayGamma_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaLarge_bin, h_invXsec_decay_photons_etaLarge, file, dir_decayGamma_invXsec, 2.*etaLarge);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaTPC_bin, h_invXsec_decay_photons_etaTPC, file, dir_decayGamma_invXsec, 2.*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaEMCal_bin, h_invXsec_decay_photons_etaEMCal, file, dir_decayGamma_invXsec, 2.*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaPHOS_bin, h_invXsec_decay_photons_etaPHOS, file, dir_decayGamma_invXsec, 2.*etaPHOS);


  //-----------------------------
  file.Close();

  return 0;
}
