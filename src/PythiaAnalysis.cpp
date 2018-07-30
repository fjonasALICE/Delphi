#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "PythiaAnalysisHelper.h"
#include "PythiaAnalysis.h"
#include "TMath.h"

#include "fastjet/ClusterSequence.hh"

using std::cout;
using namespace Pythia8;
using fastjet::PseudoJet;
using fastjet::JetDefinition;
using fastjet::ClusterSequence;
using fastjet::antikt_algorithm;

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

  if (!strcmp(argv[2],"MBVeto")){ // prevent double sampling of MB events and e.g. JJ events
    MB_veto = true;
    printf("\n Applying MBVeto, i.e. generating MB events with restriction pthat < pTHatBin[0] = %f\n(in order to prevent double sampling with events generated in pthatbins)\n", pTHatBin[0]); 
  }

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

  //----------------------------------------------------------------------------------------------------
  // check underlying born kt to see, e.g., if HardQCD cross section does not blow up
  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", 1000, ptMin, ptMax);
  // store the weight sum for proper normalization afterwards
  TH1D *h_weightSum = new TH1D("h_weightSum","weightSum = number of events for pythia standalone", 1, 0., 1.);
  
  // TH2D electron_pt vs electron_topMotherID
  TH2D *h2_electron_pt_topMotherID = new TH2D("h2_electron_pt_topMotherID","electron_pt_topMotherID (EMCal acceptance |#eta| < 0.66)",17,0,17,ptBins,ptMin,ptMax);
  h2_electron_pt_topMotherID->SetCanExtend(TH1::kXaxis);
  for(int i = 0; i < 17; i++)
    h2_electron_pt_topMotherID->GetXaxis()->SetBinLabel(i+1, electronMotherName[i]);

  // charged jets
  TH1D *h_chJets_pt_etaTPC = new TH1D("h_chjet_pt_etaTPC",Form("charged jet pt in |#eta| < (0.9-R), R=%f",jetRadius), ptBins, ptMin, ptMax);
  TH1D *h_chJets_pt_leading_etaTPC = new TH1D("h_chjet_pt_leading_etaTPC",Form("leading charged jet pt in |#eta| < (0.9-R), R=%f",jetRadius), ptBins, ptMin, ptMax);
  TH1D *h_dPhiJetGamma   = new TH1D("h_dPhiJetGamma"   ,"#Delta #phi_{J#gamma}", 136, -0.1, 3.3);
  TH1D *h_xJetGamma      = new TH1D("h_xJetGamma"      ,"x_{J#gamma} = p_{T}^{Jet} / p_{T}^{#gamma}", 150, 0., 3.);
  TH1D *h_chJetTrackMult = new TH1D("h_chJetTrackMult" ,"charged track multiplicity in jet", 50, 0.5, 50.5);
  TH1D *h_xObs_pGoing    = new TH1D("h_xObs_pGoing"    ,"xObs_pGoing", 1000, 0., 0.1);
  TH1D *h_xObs_PbGoing   = new TH1D("h_xObs_PbGoing"   ,"xObs_PbGoing", 1000, 0., 0.1);

  TH1D *h_isoCone_track_dPhi = new TH1D("h_isoCone_track_dPhi","#Delta #phi between photon and iso track", 100, -1., 1.);
  TH1D *h_isoCone_track_dEta = new TH1D("h_isoCone_track_dEta","#Delta #eta between photon and iso track", 100, -1., 1.);

  TH1D *h_isoPt           = new TH1D("h_isoPt"          ,"sum of pt in iso cone", ptBins, ptMin, ptMax);
  TH1D *h_isoPt_corrected = new TH1D("h_isoPt_corrected","sum of pt in iso cone minus UE", ptBins, ptMin, ptMax);

  TH1D *h_xBjorken_1 = new TH1D("h_xBjorken_1" ," Bjorken x from pythia's function x1()", 1000, 0., 0.1);
  TH1D *h_xBjorken_2 = new TH1D("h_xBjorken_2" ," Bjorken x from pythia's function x2()", 1000, 0., 0.1);
  
  TH1D *h_xSecTriggerGamma = new TH1D("h_xSecTriggerGamma","accumulated cross section of trigger photons", 1, -0.5, 0.5);
    
  // all pions without secondary correction
  TH1D *h_pi0_yDefault = new TH1D("h_pi0_yDefault","#pi^{0} in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaLarge = new TH1D("h_pi0_etaLarge","#pi^{0} in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaTPC   = new TH1D("h_pi0_etaTPC"  ,"#pi^{0} in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaEMCal = new TH1D("h_pi0_etaEMCal","#pi^{0} in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_pi0_etaPHOS  = new TH1D("h_pi0_etaPHOS" ,"#pi^{0} in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // primary pions (with secondary correction)
  TH1D *h_pi0primary_yDefault = new TH1D("h_pi0primary_yDefault","#pi^{0} (primary) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaLarge = new TH1D("h_pi0primary_etaLarge","#pi^{0} (primary) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaTPC   = new TH1D("h_pi0primary_etaTPC"  ,"#pi^{0} (primary) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaEMCal = new TH1D("h_pi0primary_etaEMCal","#pi^{0} (primary) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_pi0primary_etaPHOS  = new TH1D("h_pi0primary_etaPHOS" ,"#pi^{0} (primary) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // eta meson
  TH1D *h_eta_yDefault = new TH1D("h_eta_yDefault","#eta in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaLarge = new TH1D("h_eta_etaLarge","#eta in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaTPC   = new TH1D("h_eta_etaTPC"  ,"#eta in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaEMCal = new TH1D("h_eta_etaEMCal","#eta in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_eta_etaPHOS  = new TH1D("h_eta_etaPHOS" ,"#eta in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // eta prime meson
  TH1D *h_etaprime_yDefault = new TH1D("h_etaprime_yDefault","#eta' in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaLarge = new TH1D("h_etaprime_etaLarge","#eta' in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaTPC   = new TH1D("h_etaprime_etaTPC"  ,"#eta' in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaEMCal = new TH1D("h_etaprime_etaEMCal","#eta' in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_etaprime_etaPHOS  = new TH1D("h_etaprime_etaPHOS" ,"#eta' in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // omega meson
  TH1D *h_omega_yDefault = new TH1D("h_omega_yDefault","#omega in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaLarge = new TH1D("h_omega_etaLarge","#omega in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaTPC   = new TH1D("h_omega_etaTPC"  ,"#omega in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaEMCal = new TH1D("h_omega_etaEMCal","#omega in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_omega_etaPHOS  = new TH1D("h_omega_etaPHOS" ,"#omega in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // direct photons (consider only direct photons)
  TH1D *h_direct_photons_yDefault = new TH1D("h_direct_photons_yDefault","direct photons in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaLarge = new TH1D("h_direct_photons_etaLarge","direct photons in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaTPC   = new TH1D("h_direct_photons_etaTPC"  ,"direct photons in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaEMCal = new TH1D("h_direct_photons_etaEMCal","direct photons in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_direct_photons_etaPHOS  = new TH1D("h_direct_photons_etaPHOS" ,"direct photons in |#eta| < 0.12", ptBins, ptMin, ptMax);
  
  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_shower_photons_yDefault = new TH1D("h_shower_photons_yDefault","shower photons (q -> q #gamma) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaLarge = new TH1D("h_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaTPC   = new TH1D("h_shower_photons_etaTPC"  ,"shower photons (q -> q #gamma) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaEMCal = new TH1D("h_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_etaPHOS  = new TH1D("h_shower_photons_etaPHOS" ,"shower photons (q -> q #gamma) in |#eta| < 0.12", ptBins, ptMin, ptMax);
  
  // photons from ME (aka prompt)
  TH1D *h_222_photons_yDefault = new TH1D("h_222_photons_yDefault","photons from ME (aka prompt) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaLarge = new TH1D("h_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaTPC   = new TH1D("h_222_photons_etaTPC"  ,"photons from ME (aka prompt) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaEMCal = new TH1D("h_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_222_photons_etaPHOS  = new TH1D("h_222_photons_etaPHOS" ,"photons from ME (aka prompt) in |#eta| < 0.12", ptBins, ptMin, ptMax);
  
  // decay photons
  TH1D *h_decay_photons_yDefault = new TH1D("h_decay_photons_yDefault","decay photons in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaLarge = new TH1D("h_decay_photons_etaLarge","decay photons in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaTPC   = new TH1D("h_decay_photons_etaTPC"  ,"decay photons in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaEMCal = new TH1D("h_decay_photons_etaEMCal","decay photons in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaPHOS  = new TH1D("h_decay_photons_etaPHOS" ,"decay photons in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // isolated photons (considers only direct photons)
  TH1D *h_iso_charged2GeV_R03_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R03_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R04_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R04_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged2GeV_R05_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged2GeV_R05_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R03_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R03_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R04_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R04_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_charged3GeV_R05_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_charged3GeV_R05_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R03_photons_etaTPC      = new TH1D("h_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R03_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R04_photons_etaTPC      = new TH1D("h_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R04_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full2GeV_R05_photons_etaTPC      = new TH1D("h_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full2GeV_R05_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R03_photons_etaTPC      = new TH1D("h_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R03_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R04_photons_etaTPC      = new TH1D("h_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R04_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_iso_full3GeV_R05_photons_etaTPC      = new TH1D("h_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_iso_full3GeV_R05_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  
  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------
  // all pions without secondary correction
  TH1D *h_invXsec_pi0_yDefault = new TH1D("h_invXsec_pi0_yDefault","#pi^{0} in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaLarge = new TH1D("h_invXsec_pi0_etaLarge","#pi^{0} in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaTPC   = new TH1D("h_invXsec_pi0_etaTPC"  ,"#pi^{0} in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaEMCal = new TH1D("h_invXsec_pi0_etaEMCal","#pi^{0} in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0_etaPHOS  = new TH1D("h_invXsec_pi0_etaPHOS" ,"#pi^{0} in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // primary pions (with secondary correction)
  TH1D *h_invXsec_pi0primary_yDefault = new TH1D("h_invXsec_pi0primary_yDefault","#pi^{0} (primary) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaLarge = new TH1D("h_invXsec_pi0primary_etaLarge","#pi^{0} (primary) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaTPC   = new TH1D("h_invXsec_pi0primary_etaTPC"  ,"#pi^{0} (primary) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaEMCal = new TH1D("h_invXsec_pi0primary_etaEMCal","#pi^{0} (primary) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_pi0primary_etaPHOS  = new TH1D("h_invXsec_pi0primary_etaPHOS" ,"#pi^{0} (primary) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // eta meson
  TH1D *h_invXsec_eta_yDefault = new TH1D("h_invXsec_eta_yDefault","#eta in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaLarge = new TH1D("h_invXsec_eta_etaLarge","#eta in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaTPC   = new TH1D("h_invXsec_eta_etaTPC"  ,"#eta in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaEMCal = new TH1D("h_invXsec_eta_etaEMCal","#eta in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_eta_etaPHOS  = new TH1D("h_invXsec_eta_etaPHOS" ,"#eta in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // eta prime meson
  TH1D *h_invXsec_etaprime_yDefault = new TH1D("h_invXsec_etaprime_yDefault","#eta' in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaLarge = new TH1D("h_invXsec_etaprime_etaLarge","#eta' in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaTPC   = new TH1D("h_invXsec_etaprime_etaTPC"  ,"#eta' in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaEMCal = new TH1D("h_invXsec_etaprime_etaEMCal","#eta' in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_etaprime_etaPHOS  = new TH1D("h_invXsec_etaprime_etaPHOS" ,"#eta' in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // omega meson
  TH1D *h_invXsec_omega_yDefault = new TH1D("h_invXsec_omega_yDefault","#omega in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaLarge = new TH1D("h_invXsec_omega_etaLarge","#omega in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaTPC   = new TH1D("h_invXsec_omega_etaTPC"  ,"#omega in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaEMCal = new TH1D("h_invXsec_omega_etaEMCal","#omega in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_omega_etaPHOS  = new TH1D("h_invXsec_omega_etaPHOS" ,"#omega in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // direct photons (consider only direct photons)
  TH1D *h_invXsec_direct_photons_yDefault = new TH1D("h_invXsec_direct_photons_yDefault","direct photons in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaLarge = new TH1D("h_invXsec_direct_photons_etaLarge","direct photons in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaTPC   = new TH1D("h_invXsec_direct_photons_etaTPC"  ,"direct photons in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaEMCal = new TH1D("h_invXsec_direct_photons_etaEMCal","direct photons in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_direct_photons_etaPHOS  = new TH1D("h_invXsec_direct_photons_etaPHOS" ,"direct photons in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_invXsec_shower_photons_yDefault = new TH1D("h_invXsec_shower_photons_yDefault","shower photons (q -> q #gamma) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaLarge = new TH1D("h_invXsec_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaTPC   = new TH1D("h_invXsec_shower_photons_etaTPC"  ,"shower photons (q -> q #gamma) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaEMCal = new TH1D("h_invXsec_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_shower_photons_etaPHOS  = new TH1D("h_invXsec_shower_photons_etaPHOS" ,"shower photons (q -> q #gamma) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // photons from ME (aka prompt)
  TH1D *h_invXsec_222_photons_yDefault = new TH1D("h_invXsec_222_photons_yDefault","photons from ME (aka prompt) in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaLarge = new TH1D("h_invXsec_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaTPC   = new TH1D("h_invXsec_222_photons_etaTPC"  ,"photons from ME (aka prompt) in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaEMCal = new TH1D("h_invXsec_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_222_photons_etaPHOS  = new TH1D("h_invXsec_222_photons_etaPHOS" ,"photons from ME (aka prompt) in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // decay photons
  TH1D *h_invXsec_decay_photons_yDefault = new TH1D("h_invXsec_decay_photons_yDefault","decay photons in |y| < 0.8", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaLarge = new TH1D("h_invXsec_decay_photons_etaLarge","decay photons in |#eta| < 3.0", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaTPC   = new TH1D("h_invXsec_decay_photons_etaTPC"  ,"decay photons in |#eta| < 0.9", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaEMCal = new TH1D("h_invXsec_decay_photons_etaEMCal","decay photons in |#eta| < 0.66", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_decay_photons_etaPHOS  = new TH1D("h_invXsec_decay_photons_etaPHOS" ,"decay photons in |#eta| < 0.12", ptBins, ptMin, ptMax);

  // isolated photons (considers only direct photons)
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", ptBins, ptMin, ptMax);



  vector <TH1D*> vec_pTHat_bin;
  vector <TH1D*> vec_weightSum_bin;

  // charged jets
  vector <TH1D*> vec_chJets_pt_etaTPC_bin;
  vector <TH1D*> vec_chJets_pt_leading_etaTPC_bin;
  vector <TH1D*> vec_dPhiJetGamma_bin;
  vector <TH1D*> vec_xJetGamma_bin;
  vector <TH1D*> vec_chJetTrackMult_bin;
  vector <TH1D*> vec_xObs_pGoing_bin;
  vector <TH1D*> vec_xObs_PbGoing_bin;

  vector <TH1D*> vec_isoCone_track_dPhi_bin;
  vector <TH1D*> vec_isoCone_track_dEta_bin;

  vector <TH1D*> vec_isoPt_bin;
  vector <TH1D*> vec_isoPt_corrected_bin;

  vector <TH1D*> vec_xBjorken_1_bin;
  vector <TH1D*> vec_xBjorken_2_bin;

  vector <TH1D*> vec_xSecTriggerGamma_bin;

  // organise pTHat wise histograms in vectors
  vector <TH2D*> vec_electron_pt_topMotherID_bin;

  vector <TH1D*> vec_pi0_yDefault_bin;
  vector <TH1D*> vec_pi0_etaLarge_bin;
  vector <TH1D*> vec_pi0_etaTPC_bin;
  vector <TH1D*> vec_pi0_etaEMCal_bin;
  vector <TH1D*> vec_pi0_etaPHOS_bin;

  vector <TH1D*> vec_pi0primary_yDefault_bin;
  vector <TH1D*> vec_pi0primary_etaLarge_bin;
  vector <TH1D*> vec_pi0primary_etaTPC_bin;
  vector <TH1D*> vec_pi0primary_etaEMCal_bin;
  vector <TH1D*> vec_pi0primary_etaPHOS_bin;

  vector <TH1D*> vec_eta_yDefault_bin;
  vector <TH1D*> vec_eta_etaLarge_bin;
  vector <TH1D*> vec_eta_etaTPC_bin;
  vector <TH1D*> vec_eta_etaEMCal_bin;
  vector <TH1D*> vec_eta_etaPHOS_bin;

  vector <TH1D*> vec_etaprime_yDefault_bin;
  vector <TH1D*> vec_etaprime_etaLarge_bin;
  vector <TH1D*> vec_etaprime_etaTPC_bin;
  vector <TH1D*> vec_etaprime_etaEMCal_bin;
  vector <TH1D*> vec_etaprime_etaPHOS_bin;

  vector <TH1D*> vec_omega_yDefault_bin;
  vector <TH1D*> vec_omega_etaLarge_bin;
  vector <TH1D*> vec_omega_etaTPC_bin;
  vector <TH1D*> vec_omega_etaEMCal_bin;
  vector <TH1D*> vec_omega_etaPHOS_bin;

  vector <TH1D*> vec_direct_photons_yDefault_bin;
  vector <TH1D*> vec_direct_photons_etaLarge_bin;
  vector <TH1D*> vec_direct_photons_etaTPC_bin;
  vector <TH1D*> vec_direct_photons_etaEMCal_bin;
  vector <TH1D*> vec_direct_photons_etaPHOS_bin;

  vector <TH1D*> vec_shower_photons_yDefault_bin;
  vector <TH1D*> vec_shower_photons_etaLarge_bin;
  vector <TH1D*> vec_shower_photons_etaTPC_bin;
  vector <TH1D*> vec_shower_photons_etaEMCal_bin;
  vector <TH1D*> vec_shower_photons_etaPHOS_bin;

  vector <TH1D*> vec_222_photons_yDefault_bin;
  vector <TH1D*> vec_222_photons_etaLarge_bin;
  vector <TH1D*> vec_222_photons_etaTPC_bin;
  vector <TH1D*> vec_222_photons_etaEMCal_bin;
  vector <TH1D*> vec_222_photons_etaPHOS_bin;

  vector <TH1D*> vec_decay_photons_yDefault_bin;
  vector <TH1D*> vec_decay_photons_etaLarge_bin;
  vector <TH1D*> vec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_decay_photons_etaPHOS_bin;

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

  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------
  vector <TH1D*> vec_invXsec_pi0_yDefault_bin;
  vector <TH1D*> vec_invXsec_pi0_etaLarge_bin;
  vector <TH1D*> vec_invXsec_pi0_etaTPC_bin;
  vector <TH1D*> vec_invXsec_pi0_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_pi0_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_pi0primary_yDefault_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaLarge_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaTPC_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_pi0primary_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_eta_yDefault_bin;
  vector <TH1D*> vec_invXsec_eta_etaLarge_bin;
  vector <TH1D*> vec_invXsec_eta_etaTPC_bin;
  vector <TH1D*> vec_invXsec_eta_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_eta_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_etaprime_yDefault_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaLarge_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaTPC_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_etaprime_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_omega_yDefault_bin;
  vector <TH1D*> vec_invXsec_omega_etaLarge_bin;
  vector <TH1D*> vec_invXsec_omega_etaTPC_bin;
  vector <TH1D*> vec_invXsec_omega_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_omega_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_direct_photons_yDefault_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_direct_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_shower_photons_yDefault_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_shower_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_222_photons_yDefault_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_222_photons_etaPHOS_bin;

  vector <TH1D*> vec_invXsec_decay_photons_yDefault_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaLarge_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_invXsec_decay_photons_etaPHOS_bin;

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

  //----------------------------------------------------------------------------------------------------

  for(int i = 0; i < pTHatBins; i++){

    vec_electron_pt_topMotherID_bin.push_back( (TH2D*)h2_electron_pt_topMotherID->Clone(Form( "h2_electron_pt_topMotherID_bin_%02d", i)) );

    vec_weightSum_bin.push_back( (TH1D*)h_weightSum->Clone(Form( "h_weightSum_bin_%02d", i )) );

    // charged jets
    vec_chJets_pt_etaTPC_bin.push_back( (TH1D*)h_chJets_pt_etaTPC->Clone(Form("h_chJets_pt_etaTPC_bin_%02d",i)) );
    vec_chJets_pt_leading_etaTPC_bin.push_back( (TH1D*)h_chJets_pt_leading_etaTPC->Clone(Form("h_chJets_pt_leading_etaTPC_bin_%02d",i)) );
    vec_dPhiJetGamma_bin.push_back( (TH1D*)h_dPhiJetGamma->Clone(Form("h_dPhiJetGamma_bin_%02d",i)) );
    vec_xJetGamma_bin.push_back( (TH1D*)h_xJetGamma->Clone(Form("h_xJetGamma_bin_%02d",i)) );
    vec_chJetTrackMult_bin.push_back( (TH1D*)h_chJetTrackMult->Clone(Form("h_chJetTrackMult_bin_%02d",i)) );
    vec_xObs_pGoing_bin.push_back( (TH1D*)h_xObs_pGoing->Clone(Form("h_xObs_pGoing_bin_%02d",i)) );
    vec_xObs_PbGoing_bin.push_back( (TH1D*)h_xObs_PbGoing->Clone(Form("h_xObs_PbGoing_bin_%02d",i)) );

    vec_isoCone_track_dPhi_bin.push_back( (TH1D*)h_isoCone_track_dPhi->Clone(Form("h_isoCone_track_dPhi_bin_%02d",i)) );
    vec_isoCone_track_dEta_bin.push_back( (TH1D*)h_isoCone_track_dEta->Clone(Form("h_isoCone_track_dEta_bin_%02d",i)) );

    vec_isoPt_bin.push_back( (TH1D*)h_isoPt->Clone(Form("h_isoPt_bin_%02d",i)) );
    vec_isoPt_corrected_bin.push_back( (TH1D*)h_isoPt_corrected->Clone(Form("h_isoPt_corrected_bin_%02d",i)) );

    vec_xBjorken_1_bin.push_back( (TH1D*)h_xBjorken_1->Clone(Form("h_xBjorken_1_bin_%02d",i)) );
    vec_xBjorken_2_bin.push_back( (TH1D*)h_xBjorken_2->Clone(Form("h_xBjorken_2_bin_%02d",i)) );

    vec_xSecTriggerGamma_bin.push_back( (TH1D*)h_xSecTriggerGamma->Clone(Form( "h_xSecTriggerGamma_bin_%02d", i )) );

    // pi0 histos
    vec_pi0_yDefault_bin.push_back( (TH1D*)h_pi0_yDefault->Clone(Form("h_pi0_yDefault_bin_%02d",i)) );
    vec_pi0_etaLarge_bin.push_back( (TH1D*)h_pi0_etaLarge->Clone(Form("h_pi0_etaLarge_bin_%02d",i)) );
    vec_pi0_etaTPC_bin.push_back( (TH1D*)h_pi0_etaTPC->Clone(Form("h_pi0_etaTPC_bin_%02d",i)) );
    vec_pi0_etaEMCal_bin.push_back( (TH1D*)h_pi0_etaEMCal->Clone(Form("h_pi0_etaEMCal_bin_%02d",i)) );
    vec_pi0_etaPHOS_bin.push_back( (TH1D*)h_pi0_etaPHOS->Clone(Form("h_pi0_etaPHOS_bin_%02d",i)) );

    // primary pi0 histos
    vec_pi0primary_yDefault_bin.push_back( (TH1D*)h_pi0primary_yDefault->Clone(Form("h_pi0primary_yDefault_bin_%02d",i)) );
    vec_pi0primary_etaLarge_bin.push_back( (TH1D*)h_pi0primary_etaLarge->Clone(Form("h_pi0primary_etaLarge_bin_%02d",i)) );
    vec_pi0primary_etaTPC_bin.push_back( (TH1D*)h_pi0primary_etaTPC->Clone(Form("h_pi0primary_etaTPC_bin_%02d",i)) );
    vec_pi0primary_etaEMCal_bin.push_back( (TH1D*)h_pi0primary_etaEMCal->Clone(Form("h_pi0primary_etaEMCal_bin_%02d",i)) );
    vec_pi0primary_etaPHOS_bin.push_back( (TH1D*)h_pi0primary_etaPHOS->Clone(Form("h_pi0primary_etaPHOS_bin_%02d",i)) );

    // eta histos
    vec_eta_yDefault_bin.push_back( (TH1D*)h_eta_yDefault->Clone(Form("h_eta_yDefault_bin_%02d",i)) );
    vec_eta_etaLarge_bin.push_back( (TH1D*)h_eta_etaLarge->Clone(Form("h_eta_etaLarge_bin_%02d",i)) );
    vec_eta_etaTPC_bin.push_back( (TH1D*)h_eta_etaTPC->Clone(Form("h_eta_etaTPC_bin_%02d",i)) );
    vec_eta_etaEMCal_bin.push_back( (TH1D*)h_eta_etaEMCal->Clone(Form("h_eta_etaEMCal_bin_%02d",i)) );
    vec_eta_etaPHOS_bin.push_back( (TH1D*)h_eta_etaPHOS->Clone(Form("h_eta_etaPHOS_bin_%02d",i)) );

    // eta prime histos
    vec_etaprime_yDefault_bin.push_back( (TH1D*)h_etaprime_yDefault->Clone(Form("h_etaprime_yDefault_bin_%02d",i)) );
    vec_etaprime_etaLarge_bin.push_back( (TH1D*)h_etaprime_etaLarge->Clone(Form("h_etaprime_etaLarge_bin_%02d",i)) );
    vec_etaprime_etaTPC_bin.push_back( (TH1D*)h_etaprime_etaTPC->Clone(Form("h_etaprime_etaTPC_bin_%02d",i)) );
    vec_etaprime_etaEMCal_bin.push_back( (TH1D*)h_etaprime_etaEMCal->Clone(Form("h_etaprime_etaEMCal_bin_%02d",i)) );
    vec_etaprime_etaPHOS_bin.push_back( (TH1D*)h_etaprime_etaPHOS->Clone(Form("h_etaprime_etaPHOS_bin_%02d",i)) );

    // omega histos
    vec_omega_yDefault_bin.push_back( (TH1D*)h_omega_yDefault->Clone(Form("h_omega_yDefault_bin_%02d",i)) );
    vec_omega_etaLarge_bin.push_back( (TH1D*)h_omega_etaLarge->Clone(Form("h_omega_etaLarge_bin_%02d",i)) );
    vec_omega_etaTPC_bin.push_back( (TH1D*)h_omega_etaTPC->Clone(Form("h_omega_etaTPC_bin_%02d",i)) );
    vec_omega_etaEMCal_bin.push_back( (TH1D*)h_omega_etaEMCal->Clone(Form("h_omega_etaEMCal_bin_%02d",i)) );
    vec_omega_etaPHOS_bin.push_back( (TH1D*)h_omega_etaPHOS->Clone(Form("h_omega_etaPHOS_bin_%02d",i)) );

    // direct photon histos
    vec_direct_photons_yDefault_bin.push_back( (TH1D*)h_direct_photons_yDefault->Clone(Form("h_direct_photons_yDefault_bin_%02d",i)) );
    vec_direct_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_direct_photons_etaLarge_bin_%02d",i)) );
    vec_direct_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_direct_photons_etaTPC_bin_%02d",i)) );
    vec_direct_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_direct_photons_etaEMCal_bin_%02d",i)) );
    vec_direct_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_direct_photons_etaPHOS_bin_%02d",i)) );

    // discriminating shower and prompt photons
    vec_shower_photons_yDefault_bin.push_back( (TH1D*)h_direct_photons_yDefault->Clone(Form("h_shower_photons_yDefault_bin_%02d",i)) );
    vec_shower_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_shower_photons_etaLarge_bin_%02d",i)) );
    vec_shower_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_shower_photons_etaTPC_bin_%02d",i)) );
    vec_shower_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_shower_photons_etaEMCal_bin_%02d",i)) );
    vec_shower_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_shower_photons_etaPHOS_bin_%02d",i)) );

    vec_222_photons_yDefault_bin.push_back( (TH1D*)h_direct_photons_yDefault->Clone(Form("h_222_photons_yDefault_bin_%02d",i)) );
    vec_222_photons_etaLarge_bin.push_back( (TH1D*)h_direct_photons_etaLarge->Clone(Form("h_222_photons_etaLarge_bin_%02d",i)) );
    vec_222_photons_etaTPC_bin.push_back( (TH1D*)h_direct_photons_etaTPC->Clone(Form("h_222_photons_etaTPC_bin_%02d",i)) );
    vec_222_photons_etaEMCal_bin.push_back( (TH1D*)h_direct_photons_etaEMCal->Clone(Form("h_222_photons_etaEMCal_bin_%02d",i)) );
    vec_222_photons_etaPHOS_bin.push_back( (TH1D*)h_direct_photons_etaPHOS->Clone(Form("h_222_photons_etaPHOS_bin_%02d",i)) );

    // decay photon histos
    vec_decay_photons_yDefault_bin.push_back( (TH1D*)h_decay_photons_yDefault->Clone(Form("h_decay_photons_yDefault_bin_%02d",i)) );
    vec_decay_photons_etaLarge_bin.push_back( (TH1D*)h_decay_photons_etaLarge->Clone(Form("h_decay_photons_etaLarge_bin_%02d",i)) );
    vec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(Form("h_decay_photons_etaTPC_bin_%02d",i)) );
    vec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(Form("h_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(Form("h_decay_photons_etaPHOS_bin_%02d",i)) );

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


    //----------------------------------------------------------------------------------------------------
    // do the same jazz for invariant cross section histos ------------------------------------------
    //------------------------------------------------------------------------------------------
    // pi0 histos
    vec_invXsec_pi0_yDefault_bin.push_back( (TH1D*)h_invXsec_pi0_yDefault->Clone(Form("h_invXsec_pi0_yDefault_bin_%02d",i)) );
    vec_invXsec_pi0_etaLarge_bin.push_back( (TH1D*)h_invXsec_pi0_etaLarge->Clone(Form("h_invXsec_pi0_etaLarge_bin_%02d",i)) );
    vec_invXsec_pi0_etaTPC_bin.push_back( (TH1D*)h_invXsec_pi0_etaTPC->Clone(Form("h_invXsec_pi0_etaTPC_bin_%02d",i)) );
    vec_invXsec_pi0_etaEMCal_bin.push_back( (TH1D*)h_invXsec_pi0_etaEMCal->Clone(Form("h_invXsec_pi0_etaEMCal_bin_%02d",i)) );
    vec_invXsec_pi0_etaPHOS_bin.push_back( (TH1D*)h_invXsec_pi0_etaPHOS->Clone(Form("h_invXsec_pi0_etaPHOS_bin_%02d",i)) );

    // primary pi0 histos
    vec_invXsec_pi0primary_yDefault_bin.push_back( (TH1D*)h_invXsec_pi0primary_yDefault->Clone(Form("h_invXsec_pi0primary_yDefault_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaLarge_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaLarge->Clone(Form("h_invXsec_pi0primary_etaLarge_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaTPC_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaTPC->Clone(Form("h_invXsec_pi0primary_etaTPC_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaEMCal_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaEMCal->Clone(Form("h_invXsec_pi0primary_etaEMCal_bin_%02d",i)) );
    vec_invXsec_pi0primary_etaPHOS_bin.push_back( (TH1D*)h_invXsec_pi0primary_etaPHOS->Clone(Form("h_invXsec_pi0primary_etaPHOS_bin_%02d",i)) );

    // eta histos
    vec_invXsec_eta_yDefault_bin.push_back( (TH1D*)h_invXsec_eta_yDefault->Clone(Form("h_invXsec_eta_yDefault_bin_%02d",i)) );
    vec_invXsec_eta_etaLarge_bin.push_back( (TH1D*)h_invXsec_eta_etaLarge->Clone(Form("h_invXsec_eta_etaLarge_bin_%02d",i)) );
    vec_invXsec_eta_etaTPC_bin.push_back( (TH1D*)h_invXsec_eta_etaTPC->Clone(Form("h_invXsec_eta_etaTPC_bin_%02d",i)) );
    vec_invXsec_eta_etaEMCal_bin.push_back( (TH1D*)h_invXsec_eta_etaEMCal->Clone(Form("h_invXsec_eta_etaEMCal_bin_%02d",i)) );
    vec_invXsec_eta_etaPHOS_bin.push_back( (TH1D*)h_invXsec_eta_etaPHOS->Clone(Form("h_invXsec_eta_etaPHOS_bin_%02d",i)) );

    // eta prime histos
    vec_invXsec_etaprime_yDefault_bin.push_back( (TH1D*)h_invXsec_etaprime_yDefault->Clone(Form("h_invXsec_etaprime_yDefault_bin_%02d",i)) );
    vec_invXsec_etaprime_etaLarge_bin.push_back( (TH1D*)h_invXsec_etaprime_etaLarge->Clone(Form("h_invXsec_etaprime_etaLarge_bin_%02d",i)) );
    vec_invXsec_etaprime_etaTPC_bin.push_back( (TH1D*)h_invXsec_etaprime_etaTPC->Clone(Form("h_invXsec_etaprime_etaTPC_bin_%02d",i)) );
    vec_invXsec_etaprime_etaEMCal_bin.push_back( (TH1D*)h_invXsec_etaprime_etaEMCal->Clone(Form("h_invXsec_etaprime_etaEMCal_bin_%02d",i)) );
    vec_invXsec_etaprime_etaPHOS_bin.push_back( (TH1D*)h_invXsec_etaprime_etaPHOS->Clone(Form("h_invXsec_etaprime_etaPHOS_bin_%02d",i)) );

    // omega histos
    vec_invXsec_omega_yDefault_bin.push_back( (TH1D*)h_invXsec_omega_yDefault->Clone(Form("h_invXsec_omega_yDefault_bin_%02d",i)) );
    vec_invXsec_omega_etaLarge_bin.push_back( (TH1D*)h_invXsec_omega_etaLarge->Clone(Form("h_invXsec_omega_etaLarge_bin_%02d",i)) );
    vec_invXsec_omega_etaTPC_bin.push_back( (TH1D*)h_invXsec_omega_etaTPC->Clone(Form("h_invXsec_omega_etaTPC_bin_%02d",i)) );
    vec_invXsec_omega_etaEMCal_bin.push_back( (TH1D*)h_invXsec_omega_etaEMCal->Clone(Form("h_invXsec_omega_etaEMCal_bin_%02d",i)) );
    vec_invXsec_omega_etaPHOS_bin.push_back( (TH1D*)h_invXsec_omega_etaPHOS->Clone(Form("h_invXsec_omega_etaPHOS_bin_%02d",i)) );

    // direct photon histos
    vec_invXsec_direct_photons_yDefault_bin.push_back( (TH1D*)h_invXsec_direct_photons_yDefault->Clone(Form("h_invXsec_direct_photons_yDefault_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_direct_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_direct_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_direct_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_direct_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_direct_photons_etaPHOS_bin_%02d",i)) );

    // discriminating shower and prompt photons
    vec_invXsec_shower_photons_yDefault_bin.push_back( (TH1D*)h_invXsec_direct_photons_yDefault->Clone(Form("h_invXsec_shower_photons_yDefault_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_shower_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_shower_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_shower_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_shower_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_shower_photons_etaPHOS_bin_%02d",i)) );

    vec_invXsec_222_photons_yDefault_bin.push_back( (TH1D*)h_invXsec_direct_photons_yDefault->Clone(Form("h_invXsec_222_photons_yDefault_bin_%02d",i)) );
    vec_invXsec_222_photons_etaLarge_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaLarge->Clone(Form("h_invXsec_222_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_222_photons_etaTPC_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaTPC->Clone(Form("h_invXsec_222_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_222_photons_etaEMCal_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaEMCal->Clone(Form("h_invXsec_222_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_222_photons_etaPHOS_bin.push_back( (TH1D*)h_invXsec_direct_photons_etaPHOS->Clone(Form("h_invXsec_222_photons_etaPHOS_bin_%02d",i)) );

    // decay photon histos
    vec_invXsec_decay_photons_yDefault_bin.push_back( (TH1D*)h_decay_photons_yDefault->Clone(Form("h_invXsec_decay_photons_yDefault_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaLarge_bin.push_back( (TH1D*)h_decay_photons_etaLarge->Clone(Form("h_invXsec_decay_photons_etaLarge_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(Form("h_invXsec_decay_photons_etaTPC_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(Form("h_invXsec_decay_photons_etaEMCal_bin_%02d",i)) );
    vec_invXsec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(Form("h_invXsec_decay_photons_etaPHOS_bin_%02d",i)) );


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



      
      //------------------------------------------------------------------------------------------
      //----- jets + photon correlation ----------------------------------------------------------
      //------------------------------------------------------------------------------------------
      // reset jets
      std::vector<PseudoJet> vJets;
      std::vector<PseudoJet> vPseudo;
      ClusterSequence *cs = 0;
      if(useChargedJetsGammaCorrelations){
	for (int i = 5; i < p.event.size(); i++) {
	  if (p.event[i].isFinal() && p.event[i].isCharged()) {
	    if (TMath::Abs(p.event[i].eta()) < etaTPC){
	      vPseudo.push_back(PseudoJet(p.event[i].px(),p.event[i].py(),p.event[i].pz(),p.event[i].e()));
	    }
	  }
	}
	cs = new ClusterSequence(vPseudo, jetDef_miguel);
	vJets = sorted_by_pt(cs->inclusive_jets(2.)); // argument: jetminpt
	// loop over charged jets
	//----------------------------------------------------------------------
	if (vJets.size() != 0) {
	  for(unsigned int j = 0; j < vJets.size(); j++){
	    if(TMath::Abs(vJets.at(j).eta()) > (etaTPC-jetRadius)) continue;
	    vec_chJets_pt_etaTPC_bin.at(iBin)->Fill(vJets.at(j).pt());
	    if(j == 0)
	      vec_chJets_pt_leading_etaTPC_bin.at(iBin)->Fill(vJets.at(j).pt());
	  }
	}	
      }

      photonPtMax  = -1.;
      photonPtTemp = -1.;
      iPhoton = -1;
      // search for hardest photon in this event
      //----------------------------------------------------------------------
      for (int i = 5; i < p.event.size(); i++) {
        if (p.event[i].id() == 22 && p.event[i].isFinal() && // final photon
            p.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
            TMath::Abs(p.event[i].eta()) < (etaTPC-jetRadius)){// in maximal TPC-minus-iso-cone-radius acceptance
	  
      	  // find photonPtMax
      	  photonPtTemp = p.event[i].pT();
      	  if (photonPtTemp > photonPtMax) {
      	    photonPtMax = photonPtTemp;
      	    iPhoton = i; // remember index of hardest photon
      	  }
      	}
      }

      // loop over all direct photons
      //----------------------------------------------------------------------
      for (int i = 5; i < p.event.size(); i++) {
	bool isPhotonIsolated;
	if (p.event[i].id() == 22 && p.event[i].isFinal() && // final photon
	    p.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
	    TMath::Abs(p.event[i].eta()) < etaTPC-jetRadius){       // in maximal TPC-minus-iso-cone-radius acceptance

          // photon as pseudojet for analysis
          PseudoJet photonJet(p.event[i].px(), p.event[i].py(), p.event[i].pz(), p.event[i].e());
	  if(photonJet.pt() < 15.) continue;
	  if(photonJet.pt() > 30.) continue;
	  // calculate ue pt density for a given photon i NOT IMPLEMENTED IN THE MOMENT
	  double UEPtDensity = 0.;
	  // printf("UEPtDensity(p.event, i) = %f\n",UEPtDensity);
	  // check isolation
	  isPhotonIsolated = pyHelp.IsPhotonIsolated(p.event, i, etaTPC, isoConeRadius, isoPtMax, UEPtDensity, vec_isoCone_track_dPhi_bin.at(iBin), vec_isoCone_track_dEta_bin.at(iBin), vec_isoPt_bin.at(iBin), vec_isoPt_corrected_bin.at(iBin));

	  // Fill histograms
	  //----------------------------------------------------------------------
	  // vec_directphoton_pt_bin.at(iBin)->Fill(p.event[i].pT());
	  // if(i==iPhoton) vec_directphoton_pt_leading_bin.at(iBin)->Fill(p.event[i].pT());

	  // if(isPhotonIsolated){
	  //   vec_isodirectphoton_pt_bin.at(iBin)->Fill(p.event[i].pT());
	  //   if(i==iPhoton) vec_isodirectphoton_pt_leading_bin.at(iBin)->Fill(p.event[i].pT());
	  // }

	  if(vJets.size() > 0 && isPhotonIsolated)
	    for(unsigned int iJet = 0; iJet < vJets.size(); iJet++){
	      bool isJetSeparated = ( TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))) > TMath::Pi()/2. );
	      if(vJets.at(iJet).pt() < 10.) break; // vJets are sorted by pt, break is ok
	      // gamma-jet correlation	 
	      vec_dPhiJetGamma_bin.at(iBin)->Fill(TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))));
	      if(!isJetSeparated) continue;
	      // x_Jet-gamma
	      vec_xSecTriggerGamma_bin.at(iBin)->Fill(0.);
	      vector<PseudoJet> vec_jetConst = vJets.at(iJet).constituents();
	      vec_xJetGamma_bin.at(iBin)->Fill(vJets.at(iJet).pt()/photonJet.pt());
	      // charged particle multiplicity in jets
	      vec_chJetTrackMult_bin.at(iBin)->Fill(vec_jetConst.size());
	      // x_obs p-going direction
	      vec_xObs_pGoing_bin.at(iBin)->Fill(pyHelp.XObs_pGoing(vJets.at(iJet), photonJet, p.info.eB()));
	      // x_obs Pb-going direction
	      vec_xObs_PbGoing_bin.at(iBin)->Fill(pyHelp.XObs_PbGoing(vJets.at(iJet), photonJet, p.info.eB()));
	      // real Bjorken x
	      vec_xBjorken_1_bin.at(iBin)->Fill(p.info.x1());
	      vec_xBjorken_2_bin.at(iBin)->Fill(p.info.x2());
	    }

	  // print scales of event
	  // printf("---------------------------------------\n");
	  // printf("p.info.pTHat()   = %f\n", p.info.pTHat());
	  // printf("p.info.QFac()    = %f\n", p.info.QFac());
	  // printf("p.info.QRen()    = %f\n", p.info.QRen());
	  // printf("p.info.scalup()  = %f\n", p.info.scalup());
	  // printf("\n");
	  
	  
	} // if direct photon in acceptance
      } // particle for-loop
      //------------------------------------------------------------------------------------------
      //----- END OF jets + photon correlation ---------------------------------------------------
      //------------------------------------------------------------------------------------------
      delete cs;     

      //------------------------------------------------------------------------------------------
      pyHelp.Fill_TH2_Electron_TopMotherID(p.event, etaEMCal, vec_electron_pt_topMotherID_bin.at(iBin));

      pyHelp.Fill_Pi0_Pt(p.event, yDefault, true, vec_pi0_yDefault_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaLarge, false, vec_pi0_etaLarge_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaTPC, false, vec_pi0_etaTPC_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaEMCal, false, vec_pi0_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Pi0_Pt(p.event, etaPHOS, false, vec_pi0_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Pi0Primary_Pt(p.event, yDefault, true, vec_pi0primary_yDefault_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaLarge, false, vec_pi0primary_etaLarge_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaTPC, false, vec_pi0primary_etaTPC_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaEMCal, false, vec_pi0primary_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Pi0Primary_Pt(p.event, etaPHOS, false, vec_pi0primary_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Eta_Pt(p.event, yDefault, true, vec_eta_yDefault_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaLarge, false, vec_eta_etaLarge_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaTPC, false, vec_eta_etaTPC_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaEMCal, false, vec_eta_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Eta_Pt(p.event, etaPHOS, false, vec_eta_etaPHOS_bin.at(iBin));

      pyHelp.Fill_EtaPrime_Pt(p.event, yDefault, true, vec_etaprime_yDefault_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaLarge, false, vec_etaprime_etaLarge_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaTPC, false, vec_etaprime_etaTPC_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaEMCal, false, vec_etaprime_etaEMCal_bin.at(iBin));
      pyHelp.Fill_EtaPrime_Pt(p.event, etaPHOS, false, vec_etaprime_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Omega_Pt(p.event, yDefault, true, vec_omega_yDefault_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaLarge, false, vec_omega_etaLarge_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaTPC, false, vec_omega_etaTPC_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaEMCal, false, vec_omega_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Omega_Pt(p.event, etaPHOS, false, vec_omega_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Direct_Photon_Pt(p.event, yDefault, vec_direct_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaLarge, vec_direct_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaTPC, vec_direct_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaEMCal, vec_direct_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Direct_Photon_Pt(p.event, etaPHOS, vec_direct_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_Shower_Photon_Pt(p.event, yDefault, vec_shower_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaLarge, vec_shower_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaTPC, vec_shower_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaEMCal, vec_shower_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Shower_Photon_Pt(p.event, etaPHOS, vec_shower_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_222_Photon_Pt(p.event, yDefault, vec_222_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaLarge, vec_222_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaTPC, vec_222_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaEMCal, vec_222_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_222_Photon_Pt(p.event, etaPHOS, vec_222_photons_etaPHOS_bin.at(iBin));
      
      pyHelp.Fill_Decay_Photon_Pt(p.event, yDefault, vec_decay_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaLarge, vec_decay_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaTPC, vec_decay_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaEMCal, vec_decay_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Decay_Photon_Pt(p.event, etaPHOS, vec_decay_photons_etaPHOS_bin.at(iBin));

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

      //----------------------------------------------------------------------------------------------------
      // do the same jazz for invariant cross section histos ------------------------------------------
      //------------------------------------------------------------------------------------------

      pyHelp.Fill_invXsec_Pi0_Pt(p.event, yDefault, true, vec_invXsec_pi0_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaLarge, false, vec_invXsec_pi0_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaTPC, false, vec_invXsec_pi0_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaEMCal, false, vec_invXsec_pi0_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0_Pt(p.event, etaPHOS, false, vec_invXsec_pi0_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, yDefault, true, vec_invXsec_pi0primary_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaLarge, false, vec_invXsec_pi0primary_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaTPC, false, vec_invXsec_pi0primary_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaEMCal, false, vec_invXsec_pi0primary_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Pi0Primary_Pt(p.event, etaPHOS, false, vec_invXsec_pi0primary_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Eta_Pt(p.event, yDefault, true, vec_invXsec_eta_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaLarge, false, vec_invXsec_eta_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaTPC, false, vec_invXsec_eta_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaEMCal, false, vec_invXsec_eta_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Eta_Pt(p.event, etaPHOS, false, vec_invXsec_eta_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, yDefault, true, vec_invXsec_etaprime_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaLarge, false, vec_invXsec_etaprime_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaTPC, false, vec_invXsec_etaprime_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaEMCal, false, vec_invXsec_etaprime_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_EtaPrime_Pt(p.event, etaPHOS, false, vec_invXsec_etaprime_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Omega_Pt(p.event, yDefault, true, vec_invXsec_omega_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaLarge, false, vec_invXsec_omega_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaTPC, false, vec_invXsec_omega_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaEMCal, false, vec_invXsec_omega_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Omega_Pt(p.event, etaPHOS, false, vec_invXsec_omega_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, yDefault, vec_invXsec_direct_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaLarge, vec_invXsec_direct_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaTPC, vec_invXsec_direct_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaEMCal, vec_invXsec_direct_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Direct_Photon_Pt(p.event, etaPHOS, vec_invXsec_direct_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, yDefault, vec_invXsec_shower_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaLarge, vec_invXsec_shower_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaTPC, vec_invXsec_shower_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaEMCal, vec_invXsec_shower_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_Shower_Photon_Pt(p.event, etaPHOS, vec_invXsec_shower_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, yDefault, vec_invXsec_222_photons_yDefault_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaLarge, vec_invXsec_222_photons_etaLarge_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaTPC, vec_invXsec_222_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaEMCal, vec_invXsec_222_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_invXsec_222_Photon_Pt(p.event, etaPHOS, vec_invXsec_222_photons_etaPHOS_bin.at(iBin));

      pyHelp.Fill_invXsec_Decay_Photon_Pt(p.event, yDefault, vec_invXsec_decay_photons_yDefault_bin.at(iBin));
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

    vec_chJets_pt_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_chJets_pt_leading_etaTPC_bin.at(iBin)->Scale(sigma);
    
    vec_dPhiJetGamma_bin.at(iBin)->Scale(sigma);
    vec_xJetGamma_bin.at(iBin)->Scale(sigma);
    vec_chJetTrackMult_bin.at(iBin)->Scale(sigma);
    vec_xObs_pGoing_bin.at(iBin)->Scale(sigma);
    vec_xObs_PbGoing_bin.at(iBin)->Scale(sigma);
    vec_isoCone_track_dPhi_bin.at(iBin)->Scale(sigma);
    vec_isoCone_track_dEta_bin.at(iBin)->Scale(sigma);

    vec_isoPt_bin.at(iBin)->Scale(sigma);
    vec_isoPt_corrected_bin.at(iBin)->Scale(sigma);
  
    vec_xBjorken_1_bin.at(iBin)->Scale(sigma);
    vec_xBjorken_2_bin.at(iBin)->Scale(sigma);

    vec_xSecTriggerGamma_bin.at(iBin)->Scale(sigma);
    
    vec_electron_pt_topMotherID_bin.at(iBin)->Scale(sigma);

    vec_pi0_yDefault_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_pi0_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_pi0primary_yDefault_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_pi0primary_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_eta_yDefault_bin.at(iBin)->Scale(sigma);
    vec_eta_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_eta_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_eta_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_eta_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_etaprime_yDefault_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_etaprime_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_omega_yDefault_bin.at(iBin)->Scale(sigma);
    vec_omega_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_omega_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_omega_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_omega_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_direct_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_shower_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_shower_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_222_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_222_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_decay_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

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

    //----------------------------------------------------------------------------------------------------
    // do the same jazz for invariant cross section histos ------------------------------------------
    //------------------------------------------------------------------------------------------

    vec_invXsec_pi0_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_pi0primary_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_pi0primary_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_eta_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_eta_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_etaprime_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_etaprime_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_omega_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_omega_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_direct_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_shower_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_shower_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_222_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_222_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

    vec_invXsec_decay_photons_yDefault_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_invXsec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);


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

  TDirectory *dir_chJets = file.mkdir("chJets");
  pyHelp.Add_Histos_Scale_Write2File( vec_chJets_pt_etaTPC_bin, h_chJets_pt_etaTPC , file, dir_chJets, 2*(etaTPC-jetRadius), false);
  pyHelp.Add_Histos_Scale_Write2File( vec_chJets_pt_leading_etaTPC_bin, h_chJets_pt_leading_etaTPC , file, dir_chJets, 2*(etaTPC-jetRadius), false);

  pyHelp.Add_Histos_Scale_Write2File( vec_dPhiJetGamma_bin, h_dPhiJetGamma, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_xJetGamma_bin, h_xJetGamma, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_chJetTrackMult_bin, h_chJetTrackMult, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_xObs_pGoing_bin, h_xObs_pGoing, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_xObs_PbGoing_bin, h_xObs_PbGoing, file, dir_chJets, 1., false);
  
  pyHelp.Add_Histos_Scale_Write2File( vec_isoCone_track_dPhi_bin, h_isoCone_track_dPhi, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_isoCone_track_dEta_bin, h_isoCone_track_dEta, file, dir_chJets, 1., false);

  pyHelp.Add_Histos_Scale_Write2File( vec_isoPt_bin, h_isoPt, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_isoPt_corrected_bin, h_isoPt_corrected, file, dir_chJets, 1., false);
  
  pyHelp.Add_Histos_Scale_Write2File( vec_xBjorken_1_bin, h_xBjorken_1, file, dir_chJets, 1., false);
  pyHelp.Add_Histos_Scale_Write2File( vec_xBjorken_2_bin, h_xBjorken_2, file, dir_chJets, 1., false);

  pyHelp.Add_Histos_Scale_Write2File( vec_xSecTriggerGamma_bin, h_xSecTriggerGamma, file, dir_chJets, 1., false);

  TDirectory *dir_electron = file.mkdir("electron");
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_pt_topMotherID_bin, h2_electron_pt_topMotherID, file, dir_electron, 2*etaEMCal, false);

  TDirectory *dir_pTHat = file.mkdir("pTHat");
  pyHelp.Add_Histos_Scale_Write2File( vec_pTHat_bin, h_pTHat, file, dir_pTHat, 1., false);

  TDirectory *dir_pi0 = file.mkdir("pi0");
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_yDefault_bin, h_pi0_yDefault, file, dir_pi0, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaLarge_bin, h_pi0_etaLarge, file, dir_pi0, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaTPC_bin, h_pi0_etaTPC, file, dir_pi0, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaEMCal_bin, h_pi0_etaEMCal, file, dir_pi0, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0_etaPHOS_bin, h_pi0_etaPHOS, file, dir_pi0, 2*etaPHOS, false);

  TDirectory *dir_pi0primary = file.mkdir("pi0primary");
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_yDefault_bin, h_pi0primary_yDefault, file, dir_pi0primary, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaLarge_bin, h_pi0primary_etaLarge, file, dir_pi0primary, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaTPC_bin, h_pi0primary_etaTPC, file, dir_pi0primary, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaEMCal_bin, h_pi0primary_etaEMCal, file, dir_pi0primary, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_pi0primary_etaPHOS_bin, h_pi0primary_etaPHOS, file, dir_pi0primary, 2*etaPHOS, false);

  TDirectory *dir_eta = file.mkdir("eta");
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_yDefault_bin, h_eta_yDefault, file, dir_eta, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaLarge_bin, h_eta_etaLarge, file, dir_eta, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaTPC_bin, h_eta_etaTPC, file, dir_eta, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaEMCal_bin, h_eta_etaEMCal, file, dir_eta, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_eta_etaPHOS_bin, h_eta_etaPHOS, file, dir_eta, 2*etaPHOS, false);

  TDirectory *dir_etaprime = file.mkdir("etaprime");
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_yDefault_bin, h_etaprime_yDefault, file, dir_etaprime, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaLarge_bin, h_etaprime_etaLarge, file, dir_etaprime, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaTPC_bin, h_etaprime_etaTPC, file, dir_etaprime, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaEMCal_bin, h_etaprime_etaEMCal, file, dir_etaprime, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_etaprime_etaPHOS_bin, h_etaprime_etaPHOS, file, dir_etaprime, 2*etaPHOS, false);

  TDirectory *dir_omega = file.mkdir("omega");
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_yDefault_bin, h_omega_yDefault, file, dir_omega, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaLarge_bin, h_omega_etaLarge, file, dir_omega, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaTPC_bin, h_omega_etaTPC, file, dir_omega, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaEMCal_bin, h_omega_etaEMCal, file, dir_omega, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_omega_etaPHOS_bin, h_omega_etaPHOS, file, dir_omega, 2*etaPHOS, false);

  TDirectory *dir_gamma = file.mkdir("gamma");
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_yDefault_bin, h_direct_photons_yDefault, file, dir_gamma, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaLarge_bin, h_direct_photons_etaLarge, file, dir_gamma, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaTPC_bin, h_direct_photons_etaTPC, file, dir_gamma, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaEMCal_bin, h_direct_photons_etaEMCal, file, dir_gamma, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_direct_photons_etaPHOS_bin, h_direct_photons_etaPHOS, file, dir_gamma, 2*etaPHOS, false);

  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_yDefault_bin, h_shower_photons_yDefault, file, dir_gamma, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaLarge_bin, h_shower_photons_etaLarge, file, dir_gamma, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaTPC_bin, h_shower_photons_etaTPC, file, dir_gamma, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaEMCal_bin, h_shower_photons_etaEMCal, file, dir_gamma, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_shower_photons_etaPHOS_bin, h_shower_photons_etaPHOS, file, dir_gamma, 2*etaPHOS, false);

  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_yDefault_bin, h_222_photons_yDefault, file, dir_gamma, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaLarge_bin, h_222_photons_etaLarge, file, dir_gamma, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaTPC_bin, h_222_photons_etaTPC, file, dir_gamma, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaEMCal_bin, h_222_photons_etaEMCal, file, dir_gamma, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_222_photons_etaPHOS_bin, h_222_photons_etaPHOS, file, dir_gamma, 2*etaPHOS, false);

  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_yDefault_bin, h_decay_photons_yDefault, file, dir_gamma, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaLarge_bin, h_decay_photons_etaLarge, file, dir_gamma, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaTPC_bin, h_decay_photons_etaTPC, file, dir_gamma, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaEMCal_bin, h_decay_photons_etaEMCal, file, dir_gamma, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaPHOS_bin, h_decay_photons_etaPHOS, file, dir_gamma, 2*etaPHOS, false);

  if(applyPhotonIso){
    TDirectory *dir_isoGamma = file.mkdir("isoGamma");
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaTPC_bin, h_iso_charged2GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaTPC_bin, h_iso_charged2GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaTPC_bin, h_iso_charged2GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaEMCal_bin, h_iso_charged2GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaEMCal_bin, h_iso_charged2GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaEMCal_bin, h_iso_charged2GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R03_photons_etaPHOS_bin, h_iso_charged2GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R04_photons_etaPHOS_bin, h_iso_charged2GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged2GeV_R05_photons_etaPHOS_bin, h_iso_charged2GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);

    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaTPC_bin, h_iso_charged3GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaTPC_bin, h_iso_charged3GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaTPC_bin, h_iso_charged3GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaEMCal_bin, h_iso_charged3GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaEMCal_bin, h_iso_charged3GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaEMCal_bin, h_iso_charged3GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R03_photons_etaPHOS_bin, h_iso_charged3GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R04_photons_etaPHOS_bin, h_iso_charged3GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_charged3GeV_R05_photons_etaPHOS_bin, h_iso_charged3GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);


    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaTPC_bin, h_iso_full2GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaTPC_bin, h_iso_full2GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaTPC_bin, h_iso_full2GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaEMCal_bin, h_iso_full2GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaEMCal_bin, h_iso_full2GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaEMCal_bin, h_iso_full2GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R03_photons_etaPHOS_bin, h_iso_full2GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R04_photons_etaPHOS_bin, h_iso_full2GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full2GeV_R05_photons_etaPHOS_bin, h_iso_full2GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);

    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaTPC_bin, h_iso_full3GeV_R03_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaTPC_bin, h_iso_full3GeV_R04_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaTPC_bin, h_iso_full3GeV_R05_photons_etaTPC, file, dir_isoGamma, 2*etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaEMCal_bin, h_iso_full3GeV_R03_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaEMCal_bin, h_iso_full3GeV_R04_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaEMCal_bin, h_iso_full3GeV_R05_photons_etaEMCal, file, dir_isoGamma, 2*etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R03_photons_etaPHOS_bin, h_iso_full3GeV_R03_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R04_photons_etaPHOS_bin, h_iso_full3GeV_R04_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File( vec_iso_full3GeV_R05_photons_etaPHOS_bin, h_iso_full3GeV_R05_photons_etaPHOS, file, dir_isoGamma, 2*etaPHOS, false);
  }
  

  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------

  TDirectory *dir_pi0_invXsec = file.mkdir("pi0_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_yDefault_bin, h_invXsec_pi0_yDefault, file, dir_pi0_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaLarge_bin, h_invXsec_pi0_etaLarge, file, dir_pi0_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaTPC_bin, h_invXsec_pi0_etaTPC, file, dir_pi0_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaEMCal_bin, h_invXsec_pi0_etaEMCal, file, dir_pi0_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0_etaPHOS_bin, h_invXsec_pi0_etaPHOS, file, dir_pi0_invXsec, 2.*etaPHOS, false, true);

  TDirectory *dir_pi0primary_invXsec = file.mkdir("pi0primary_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_yDefault_bin, h_invXsec_pi0primary_yDefault, file, dir_pi0primary_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaLarge_bin, h_invXsec_pi0primary_etaLarge, file, dir_pi0primary_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaTPC_bin, h_invXsec_pi0primary_etaTPC, file, dir_pi0primary_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaEMCal_bin, h_invXsec_pi0primary_etaEMCal, file, dir_pi0primary_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_pi0primary_etaPHOS_bin, h_invXsec_pi0primary_etaPHOS, file, dir_pi0primary_invXsec, 2.*etaPHOS, false, true);

  TDirectory *dir_eta_invXsec = file.mkdir("eta_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_yDefault_bin, h_invXsec_eta_yDefault, file, dir_eta_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaLarge_bin, h_invXsec_eta_etaLarge, file, dir_eta_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaTPC_bin, h_invXsec_eta_etaTPC, file, dir_eta_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaEMCal_bin, h_invXsec_eta_etaEMCal, file, dir_eta_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_eta_etaPHOS_bin, h_invXsec_eta_etaPHOS, file, dir_eta_invXsec, 2.*etaPHOS, false, true);

  TDirectory *dir_etaprime_invXsec = file.mkdir("etaprime_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_yDefault_bin, h_invXsec_etaprime_yDefault, file, dir_etaprime_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaLarge_bin, h_invXsec_etaprime_etaLarge, file, dir_etaprime_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaTPC_bin, h_invXsec_etaprime_etaTPC, file, dir_etaprime_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaEMCal_bin, h_invXsec_etaprime_etaEMCal, file, dir_etaprime_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_etaprime_etaPHOS_bin, h_invXsec_etaprime_etaPHOS, file, dir_etaprime_invXsec, 2.*etaPHOS, false, true);

  TDirectory *dir_omega_invXsec = file.mkdir("omega_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_yDefault_bin, h_invXsec_omega_yDefault, file, dir_omega_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaLarge_bin, h_invXsec_omega_etaLarge, file, dir_omega_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaTPC_bin, h_invXsec_omega_etaTPC, file, dir_omega_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaEMCal_bin, h_invXsec_omega_etaEMCal, file, dir_omega_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_omega_etaPHOS_bin, h_invXsec_omega_etaPHOS, file, dir_omega_invXsec, 2.*etaPHOS, false, true);

  TDirectory *dir_gamma_invXsec = file.mkdir("gamma_invXsec");
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_yDefault_bin, h_invXsec_direct_photons_yDefault, file, dir_gamma_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaLarge_bin, h_invXsec_direct_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaTPC_bin, h_invXsec_direct_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaEMCal_bin, h_invXsec_direct_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_direct_photons_etaPHOS_bin, h_invXsec_direct_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS, false, true);

  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_yDefault_bin, h_invXsec_shower_photons_yDefault, file, dir_gamma_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaLarge_bin, h_invXsec_shower_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaTPC_bin, h_invXsec_shower_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaEMCal_bin, h_invXsec_shower_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_shower_photons_etaPHOS_bin, h_invXsec_shower_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS, false, true);

  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_yDefault_bin, h_invXsec_222_photons_yDefault, file, dir_gamma_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaLarge_bin, h_invXsec_222_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaTPC_bin, h_invXsec_222_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaEMCal_bin, h_invXsec_222_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_222_photons_etaPHOS_bin, h_invXsec_222_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS, false, true);

  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_yDefault_bin, h_invXsec_decay_photons_yDefault, file, dir_gamma_invXsec, 2.*yDefault, true, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaLarge_bin, h_invXsec_decay_photons_etaLarge, file, dir_gamma_invXsec, 2.*etaLarge, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaTPC_bin, h_invXsec_decay_photons_etaTPC, file, dir_gamma_invXsec, 2.*etaTPC, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaEMCal_bin, h_invXsec_decay_photons_etaEMCal, file, dir_gamma_invXsec, 2.*etaEMCal, false, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_decay_photons_etaPHOS_bin, h_invXsec_decay_photons_etaPHOS, file, dir_gamma_invXsec, 2.*etaPHOS, false, true);

  if(applyPhotonIso){
    TDirectory *dir_isoGamma_invXsec = file.mkdir("isoGamma_invXsec");
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaTPC_bin, h_invXsec_iso_charged2GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_charged2GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged2GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_charged2GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);

    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaTPC_bin, h_invXsec_iso_charged3GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_charged3GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_charged3GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_charged3GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);


    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaTPC_bin, h_invXsec_iso_full2GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_full2GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full2GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_full2GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);

    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R03_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R04_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaTPC_bin, h_invXsec_iso_full3GeV_R05_photons_etaTPC, file, dir_isoGamma_invXsec, 2.*etaTPC, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R03_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R04_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaEMCal_bin, h_invXsec_iso_full3GeV_R05_photons_etaEMCal, file, dir_isoGamma_invXsec, 2.*etaEMCal, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R03_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R03_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R04_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R04_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
    pyHelp.Add_Histos_Scale_Write2File( vec_invXsec_iso_full3GeV_R05_photons_etaPHOS_bin, h_invXsec_iso_full3GeV_R05_photons_etaPHOS, file, dir_isoGamma_invXsec, 2.*etaPHOS, false, true);
  }


  //-----------------------------
  file.Close();

  return 0;
}
