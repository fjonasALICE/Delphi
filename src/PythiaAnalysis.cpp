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
#include "TCanvas.h"
#include "TPaveText.h"

#include "fastjet/ClusterSequence.hh"


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

  string pdfA, pdfB;
  if (argc >= 10){ // choice of external PDF (using LHAPDF6)
    pdfA = argv[9];
    pdfB = argv[9];
    printf("\nUsing PDF %s for beam A\n", pdfA.c_str());
    p.readString("PDF:pSet = LHAPDF6:" + pdfA);
    if( argc >= 11)
      pdfB = argv[10];
    printf("Using PDF %s for beam B\n", pdfB.c_str());
    p.readString("PDF:pSetB = LHAPDF6:" + pdfB);
  }

  int nEvent = strtol(argv[3], NULL, 10); // number of events
  if( !strcmp(argv[2],"JJ") || !strcmp(argv[2],"GJ") || !strcmp(argv[2],"WeakBoson") )
    printf("\nGenerating %d events per pthat bin\n", nEvent);
  else
    printf("\nGenerating %d events\n", nEvent);

  pyHelp.Pass_Parameters_To_Pythia(p, argc, argv); // which energy, scales, optional master switches

  int pTHatBins = 0;
  double pTHatBin[100];
  // pthat bin definition from ALICE JJ production at 8 TeV
  const int pTHatBins_250GeV = 18;
  double pTHatBin_250GeV[pTHatBins_250GeV+1] = { 9.  , 12. , 16. , 21. , 28.,
						 36. , 45. , 57. , 70. , 85.,
						 99. , 115., 132., 150., 169.,
						 190., 212., 235 , 10000. }; 

  const int pTHatBins_100GeV = 8;
  double pTHatBin_100GeV[pTHatBins_100GeV+1]= { 9.  , 12. , 16. , 21. , 28.,
						36. , 45. , 57., 10000. }; 
  // pthat bin use (ignore for MB production)
  if(useGammaJetCorrelations){
    pTHatBins = pTHatBins_100GeV;
    std::copy(pTHatBin_100GeV,pTHatBin_100GeV+9,pTHatBin);    
  }else{
    pTHatBins = pTHatBins_250GeV;
    std::copy(pTHatBin_250GeV,pTHatBin_250GeV+19,pTHatBin);    
  }

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
  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", pyHelp.ptBins, pyHelp.ptBinArray);
  // store the weight sum for proper normalization afterwards
  TH1D *h_weightSum = new TH1D("h_weightSum","weightSum = number of events for pythia standalone", 1, 0., 1.);
  
  // TH2D electron_pt vs electron_topMotherID
  TH2D *h2_electron_pt_topMotherID = new TH2D("h2_electron_pt_topMotherID","electron_pt_topMotherID (EMCal acceptance |#eta| < 0.66)",17,0,17,pyHelp.ptBins, pyHelp.ptBinArray);
  h2_electron_pt_topMotherID->SetCanExtend(TH1::kXaxis);
  for(int i = 0; i < 17; i++)
    h2_electron_pt_topMotherID->GetXaxis()->SetBinLabel(i+1, electronMotherName[i]);

  // charged jets
  TH1D *h_chJets_pt_etaTPC = new TH1D("h_chjet_pt_etaTPC",Form("charged jet pt in |#eta| < (0.9-R), R=%f",jetRadius), pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_chJets_pt_leading_etaTPC = new TH1D("h_chjet_pt_leading_etaTPC",Form("leading charged jet pt in |#eta| < (0.9-R), R=%f",jetRadius), pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_dPhiJetGamma   = new TH1D("h_dPhiJetGamma"   ,"#Delta #phi_{J#gamma}", pyHelp.dPhiJetGamma_nBins, pyHelp.dPhiJetGamma_min, pyHelp.dPhiJetGamma_max);
  TH1D *h_xJetGamma      = new TH1D("h_xJetGamma"      ,"x_{J#gamma} = p_{T}^{Jet} / p_{T}^{#gamma}", pyHelp.dxJetGamma_nBins, pyHelp.dxJetGamma_min, pyHelp.dxJetGamma_max);
  TH1D *h_chJetTrackMult = new TH1D("h_chJetTrackMult" ,"charged track multiplicity in jet", pyHelp.chJetTrackMult_nBins, pyHelp.chJetTrackMult_min, pyHelp.chJetTrackMult_max);
  TH1D *h_xObs_pGoing    = new TH1D("h_xObs_pGoing"    ,"xObs_pGoing", pyHelp.xObs_nBins, pyHelp.xObs_min, pyHelp.xObs_max);
  TH1D *h_xObs_PbGoing   = new TH1D("h_xObs_PbGoing"   ,"xObs_PbGoing", pyHelp.xObs_nBins, pyHelp.xObs_min, pyHelp.xObs_max);

  TH1D *h_xBjorken_1 = new TH1D("h_xBjorken_1" ," Bjorken x from pythia's function x1()", pyHelp.xObs_nBins, pyHelp.xObs_min, pyHelp.xObs_max);
  TH1D *h_xBjorken_2 = new TH1D("h_xBjorken_2" ," Bjorken x from pythia's function x2()", pyHelp.xObs_nBins, pyHelp.xObs_min, pyHelp.xObs_max);
  TH1D *h_isoCone_track_dPhi = new TH1D("h_isoCone_track_dPhi","#Delta #phi between photon and iso track", pyHelp.isoCone_track_nBins, pyHelp.isoCone_track_min, pyHelp.isoCone_track_max);
  TH1D *h_isoCone_track_dEta = new TH1D("h_isoCone_track_dEta","#Delta #eta between photon and iso track", pyHelp.isoCone_track_nBins, pyHelp.isoCone_track_min, pyHelp.isoCone_track_max);

  TH1D *h_xSecTriggerGamma = new TH1D("h_xSecTriggerGamma","accumulated cross section of trigger photons", 1, -0.5, 0.5);

  TH1D *h_isoPt           = new TH1D("h_isoPt"          ,"sum of pt in iso cone", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_isoPt_corrected = new TH1D("h_isoPt_corrected","sum of pt in iso cone minus UE", pyHelp.ptBins, pyHelp.ptBinArray);
    
  // all electrons (+ positrons)
  TH1D *h_electron_yDefault = new TH1D("h_electron_yDefault","e^{#pm} in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_electron_etaLarge = new TH1D("h_electron_etaLarge","e^{#pm} in |#eta| < 3.00", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_electron_etaTPC   = new TH1D("h_electron_etaTPC"  ,"e^{#pm} in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_electron_etaEMCal = new TH1D("h_electron_etaEMCal","e^{#pm} in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_electron_etaPHOS  = new TH1D("h_electron_etaPHOS" ,"e^{#pm} in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // all pions without secondary correction
  TH1D *h_pi0_yDefault = new TH1D("h_pi0_yDefault","#pi^{0} in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0_etaLarge = new TH1D("h_pi0_etaLarge","#pi^{0} in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0_etaTPC   = new TH1D("h_pi0_etaTPC"  ,"#pi^{0} in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0_etaEMCal = new TH1D("h_pi0_etaEMCal","#pi^{0} in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0_etaPHOS  = new TH1D("h_pi0_etaPHOS" ,"#pi^{0} in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // primary pions (with secondary correction)
  TH1D *h_pi0primary_yDefault = new TH1D("h_pi0primary_yDefault","#pi^{0} (primary) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0primary_etaLarge = new TH1D("h_pi0primary_etaLarge","#pi^{0} (primary) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0primary_etaTPC   = new TH1D("h_pi0primary_etaTPC"  ,"#pi^{0} (primary) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0primary_etaEMCal = new TH1D("h_pi0primary_etaEMCal","#pi^{0} (primary) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_pi0primary_etaPHOS  = new TH1D("h_pi0primary_etaPHOS" ,"#pi^{0} (primary) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // eta meson
  TH1D *h_eta_yDefault = new TH1D("h_eta_yDefault","#eta in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_eta_etaLarge = new TH1D("h_eta_etaLarge","#eta in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_eta_etaTPC   = new TH1D("h_eta_etaTPC"  ,"#eta in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_eta_etaEMCal = new TH1D("h_eta_etaEMCal","#eta in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_eta_etaPHOS  = new TH1D("h_eta_etaPHOS" ,"#eta in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // eta prime meson
  TH1D *h_etaprime_yDefault = new TH1D("h_etaprime_yDefault","#eta' in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_etaprime_etaLarge = new TH1D("h_etaprime_etaLarge","#eta' in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_etaprime_etaTPC   = new TH1D("h_etaprime_etaTPC"  ,"#eta' in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_etaprime_etaEMCal = new TH1D("h_etaprime_etaEMCal","#eta' in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_etaprime_etaPHOS  = new TH1D("h_etaprime_etaPHOS" ,"#eta' in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // omega meson
  TH1D *h_omega_yDefault = new TH1D("h_omega_yDefault","#omega in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_omega_etaLarge = new TH1D("h_omega_etaLarge","#omega in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_omega_etaTPC   = new TH1D("h_omega_etaTPC"  ,"#omega in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_omega_etaEMCal = new TH1D("h_omega_etaEMCal","#omega in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_omega_etaPHOS  = new TH1D("h_omega_etaPHOS" ,"#omega in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // direct photons (consider only direct photons)
  TH1D *h_direct_photons_yDefault = new TH1D("h_direct_photons_yDefault","direct photons in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_direct_photons_etaLarge = new TH1D("h_direct_photons_etaLarge","direct photons in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_direct_photons_etaTPC   = new TH1D("h_direct_photons_etaTPC"  ,"direct photons in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_direct_photons_etaEMCal = new TH1D("h_direct_photons_etaEMCal","direct photons in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_direct_photons_etaPHOS  = new TH1D("h_direct_photons_etaPHOS" ,"direct photons in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);
  
  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_shower_photons_yDefault = new TH1D("h_shower_photons_yDefault","shower photons (q -> q #gamma) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_shower_photons_etaLarge = new TH1D("h_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_shower_photons_etaTPC   = new TH1D("h_shower_photons_etaTPC"  ,"shower photons (q -> q #gamma) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_shower_photons_etaEMCal = new TH1D("h_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_shower_photons_etaPHOS  = new TH1D("h_shower_photons_etaPHOS" ,"shower photons (q -> q #gamma) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);
  
  // photons from ME (aka prompt)
  TH1D *h_222_photons_yDefault = new TH1D("h_222_photons_yDefault","photons from ME (aka prompt) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_222_photons_etaLarge = new TH1D("h_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_222_photons_etaTPC   = new TH1D("h_222_photons_etaTPC"  ,"photons from ME (aka prompt) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_222_photons_etaEMCal = new TH1D("h_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_222_photons_etaPHOS  = new TH1D("h_222_photons_etaPHOS" ,"photons from ME (aka prompt) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);
  
  // decay photons
  TH1D *h_decay_photons_yDefault = new TH1D("h_decay_photons_yDefault","decay photons in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_decay_photons_etaLarge = new TH1D("h_decay_photons_etaLarge","decay photons in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_decay_photons_etaTPC   = new TH1D("h_decay_photons_etaTPC"  ,"decay photons in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_decay_photons_etaEMCal = new TH1D("h_decay_photons_etaEMCal","decay photons in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_decay_photons_etaPHOS  = new TH1D("h_decay_photons_etaPHOS" ,"decay photons in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // isolated photons (considers only direct photons)
  TH1D *h_iso_charged2GeV_R03_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R03_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged2GeV_R04_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R04_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged2GeV_R05_photons_etaTPC   = new TH1D("h_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged2GeV_R05_photons_etaPHOS  = new TH1D("h_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged3GeV_R03_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R03_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged3GeV_R04_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R04_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged3GeV_R05_photons_etaTPC   = new TH1D("h_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged3GeV_R05_photons_etaPHOS  = new TH1D("h_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full2GeV_R03_photons_etaTPC      = new TH1D("h_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R03_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R03_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full2GeV_R04_photons_etaTPC      = new TH1D("h_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R04_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R04_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full2GeV_R05_photons_etaTPC      = new TH1D("h_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R05_photons_etaEMCal    = new TH1D("h_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full2GeV_R05_photons_etaPHOS     = new TH1D("h_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full3GeV_R03_photons_etaTPC      = new TH1D("h_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R03_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R03_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full3GeV_R04_photons_etaTPC      = new TH1D("h_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R04_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R04_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full3GeV_R05_photons_etaTPC      = new TH1D("h_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R05_photons_etaEMCal    = new TH1D("h_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full3GeV_R05_photons_etaPHOS     = new TH1D("h_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  // isolation sum of differen photons
  TH1D *h_iso_charged_R03_decay_photons_etaTPC = new TH1D("h_iso_charged_R03_decay_photons_etaTPC", "decay photon iso (charged pt in R=0.3) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_decay_photons_etaEMCal = new TH1D("h_iso_charged_R03_decay_photons_etaEMCal", "decay photon iso (charged pt in R=0.3) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_decay_photons_etaPHOS = new TH1D("h_iso_charged_R03_decay_photons_etaPHOS", "decay photon iso (charged pt in R=0.3) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R04_decay_photons_etaTPC = new TH1D("h_iso_charged_R04_decay_photons_etaTPC", "decay photon iso (charged pt in R=0.4) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_decay_photons_etaEMCal = new TH1D("h_iso_charged_R04_decay_photons_etaEMCal", "decay photon iso (charged pt in R=0.4) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_decay_photons_etaPHOS = new TH1D("h_iso_charged_R04_decay_photons_etaPHOS", "decay photon iso (charged pt in R=0.4) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R05_decay_photons_etaTPC = new TH1D("h_iso_charged_R05_decay_photons_etaTPC", "decay photon iso (charged pt in R=0.5) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_decay_photons_etaEMCal = new TH1D("h_iso_charged_R05_decay_photons_etaEMCal", "decay photon iso (charged pt in R=0.5) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_decay_photons_etaPHOS = new TH1D("h_iso_charged_R05_decay_photons_etaPHOS", "decay photon iso (charged pt in R=0.5) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R03_decay_photons_etaTPC = new TH1D("h_iso_full_R03_decay_photons_etaTPC", "decay photon iso (full pt in R=0.3) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_decay_photons_etaEMCal = new TH1D("h_iso_full_R03_decay_photons_etaEMCal", "decay photon iso (full pt in R=0.3) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_decay_photons_etaPHOS = new TH1D("h_iso_full_R03_decay_photons_etaPHOS", "decay photon iso (full pt in R=0.3) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R04_decay_photons_etaTPC = new TH1D("h_iso_full_R04_decay_photons_etaTPC", "decay photon iso (full pt in R=0.4) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_decay_photons_etaEMCal = new TH1D("h_iso_full_R04_decay_photons_etaEMCal", "decay photon iso (full pt in R=0.4) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_decay_photons_etaPHOS = new TH1D("h_iso_full_R04_decay_photons_etaPHOS", "decay photon iso (full pt in R=0.4) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R05_decay_photons_etaTPC = new TH1D("h_iso_full_R05_decay_photons_etaTPC", "decay photon iso (full pt in R=0.5) decay_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_decay_photons_etaEMCal = new TH1D("h_iso_full_R05_decay_photons_etaEMCal", "decay photon iso (full pt in R=0.5) decay_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_decay_photons_etaPHOS = new TH1D("h_iso_full_R05_decay_photons_etaPHOS", "decay photon iso (full pt in R=0.5) decay_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);
 // direct
  TH1D *h_iso_charged_R03_direct_photons_etaTPC = new TH1D("h_iso_charged_R03_direct_photons_etaTPC", "direct photon iso (charged pt in R=0.3) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_direct_photons_etaEMCal = new TH1D("h_iso_charged_R03_direct_photons_etaEMCal", "direct photon iso (charged pt in R=0.3) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_direct_photons_etaPHOS = new TH1D("h_iso_charged_R03_direct_photons_etaPHOS", "direct photon iso (charged pt in R=0.3) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R04_direct_photons_etaTPC = new TH1D("h_iso_charged_R04_direct_photons_etaTPC", "direct photon iso (charged pt in R=0.4) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_direct_photons_etaEMCal = new TH1D("h_iso_charged_R04_direct_photons_etaEMCal", "direct photon iso (charged pt in R=0.4) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_direct_photons_etaPHOS = new TH1D("h_iso_charged_R04_direct_photons_etaPHOS", "direct photon iso (charged pt in R=0.4) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R05_direct_photons_etaTPC = new TH1D("h_iso_charged_R05_direct_photons_etaTPC", "direct photon iso (charged pt in R=0.5) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_direct_photons_etaEMCal = new TH1D("h_iso_charged_R05_direct_photons_etaEMCal", "direct photon iso (charged pt in R=0.5) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_direct_photons_etaPHOS = new TH1D("h_iso_charged_R05_direct_photons_etaPHOS", "direct photon iso (charged pt in R=0.5) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R03_direct_photons_etaTPC = new TH1D("h_iso_full_R03_direct_photons_etaTPC", "direct photon iso (full pt in R=0.3) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_direct_photons_etaEMCal = new TH1D("h_iso_full_R03_direct_photons_etaEMCal", "direct photon iso (full pt in R=0.3) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_direct_photons_etaPHOS = new TH1D("h_iso_full_R03_direct_photons_etaPHOS", "direct photon iso (full pt in R=0.3) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R04_direct_photons_etaTPC = new TH1D("h_iso_full_R04_direct_photons_etaTPC", "direct photon iso (full pt in R=0.4) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_direct_photons_etaEMCal = new TH1D("h_iso_full_R04_direct_photons_etaEMCal", "direct photon iso (full pt in R=0.4) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_direct_photons_etaPHOS = new TH1D("h_iso_full_R04_direct_photons_etaPHOS", "direct photon iso (full pt in R=0.4) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R05_direct_photons_etaTPC = new TH1D("h_iso_full_R05_direct_photons_etaTPC", "direct photon iso (full pt in R=0.5) direct_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_direct_photons_etaEMCal = new TH1D("h_iso_full_R05_direct_photons_etaEMCal", "direct photon iso (full pt in R=0.5) direct_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_direct_photons_etaPHOS = new TH1D("h_iso_full_R05_direct_photons_etaPHOS", "direct photon iso (full pt in R=0.5) direct_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  // all
  TH1D *h_iso_charged_R03_all_photons_etaTPC = new TH1D("h_iso_charged_R03_all_photons_etaTPC", "all photon iso (charged pt in R=0.3) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_all_photons_etaEMCal = new TH1D("h_iso_charged_R03_all_photons_etaEMCal", "all photon iso (charged pt in R=0.3) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R03_all_photons_etaPHOS = new TH1D("h_iso_charged_R03_all_photons_etaPHOS", "all photon iso (charged pt in R=0.3) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R04_all_photons_etaTPC = new TH1D("h_iso_charged_R04_all_photons_etaTPC", "all photon iso (charged pt in R=0.4) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_all_photons_etaEMCal = new TH1D("h_iso_charged_R04_all_photons_etaEMCal", "all photon iso (charged pt in R=0.4) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R04_all_photons_etaPHOS = new TH1D("h_iso_charged_R04_all_photons_etaPHOS", "all photon iso (charged pt in R=0.4) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_charged_R05_all_photons_etaTPC = new TH1D("h_iso_charged_R05_all_photons_etaTPC", "all photon iso (charged pt in R=0.5) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_all_photons_etaEMCal = new TH1D("h_iso_charged_R05_all_photons_etaEMCal", "all photon iso (charged pt in R=0.5) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_charged_R05_all_photons_etaPHOS = new TH1D("h_iso_charged_R05_all_photons_etaPHOS", "all photon iso (charged pt in R=0.5) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R03_all_photons_etaTPC = new TH1D("h_iso_full_R03_all_photons_etaTPC", "all photon iso (full pt in R=0.3) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_all_photons_etaEMCal = new TH1D("h_iso_full_R03_all_photons_etaEMCal", "all photon iso (full pt in R=0.3) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R03_all_photons_etaPHOS = new TH1D("h_iso_full_R03_all_photons_etaPHOS", "all photon iso (full pt in R=0.3) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R04_all_photons_etaTPC = new TH1D("h_iso_full_R04_all_photons_etaTPC", "all photon iso (full pt in R=0.4) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_all_photons_etaEMCal = new TH1D("h_iso_full_R04_all_photons_etaEMCal", "all photon iso (full pt in R=0.4) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R04_all_photons_etaPHOS = new TH1D("h_iso_full_R04_all_photons_etaPHOS", "all photon iso (full pt in R=0.4) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_iso_full_R05_all_photons_etaTPC = new TH1D("h_iso_full_R05_all_photons_etaTPC", "all photon iso (full pt in R=0.5) all_photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_all_photons_etaEMCal = new TH1D("h_iso_full_R05_all_photons_etaEMCal", "all photon iso (full pt in R=0.5) all_photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_iso_full_R05_all_photons_etaPHOS = new TH1D("h_iso_full_R05_all_photons_etaPHOS", "all photon iso (full pt in R=0.5) all_photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);
  
  //----------------------------------------------------------------------------------------------------
  // do the same jazz for invariant cross section histos ------------------------------------------
  //------------------------------------------------------------------------------------------
  // all pions without secondary correction
  TH1D *h_invXsec_pi0_yDefault = new TH1D("h_invXsec_pi0_yDefault","#pi^{0} in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0_etaLarge = new TH1D("h_invXsec_pi0_etaLarge","#pi^{0} in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0_etaTPC   = new TH1D("h_invXsec_pi0_etaTPC"  ,"#pi^{0} in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0_etaEMCal = new TH1D("h_invXsec_pi0_etaEMCal","#pi^{0} in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0_etaPHOS  = new TH1D("h_invXsec_pi0_etaPHOS" ,"#pi^{0} in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // primary pions (with secondary correction)
  TH1D *h_invXsec_pi0primary_yDefault = new TH1D("h_invXsec_pi0primary_yDefault","#pi^{0} (primary) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0primary_etaLarge = new TH1D("h_invXsec_pi0primary_etaLarge","#pi^{0} (primary) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0primary_etaTPC   = new TH1D("h_invXsec_pi0primary_etaTPC"  ,"#pi^{0} (primary) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0primary_etaEMCal = new TH1D("h_invXsec_pi0primary_etaEMCal","#pi^{0} (primary) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_pi0primary_etaPHOS  = new TH1D("h_invXsec_pi0primary_etaPHOS" ,"#pi^{0} (primary) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // eta meson
  TH1D *h_invXsec_eta_yDefault = new TH1D("h_invXsec_eta_yDefault","#eta in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_eta_etaLarge = new TH1D("h_invXsec_eta_etaLarge","#eta in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_eta_etaTPC   = new TH1D("h_invXsec_eta_etaTPC"  ,"#eta in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_eta_etaEMCal = new TH1D("h_invXsec_eta_etaEMCal","#eta in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_eta_etaPHOS  = new TH1D("h_invXsec_eta_etaPHOS" ,"#eta in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // eta prime meson
  TH1D *h_invXsec_etaprime_yDefault = new TH1D("h_invXsec_etaprime_yDefault","#eta' in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_etaprime_etaLarge = new TH1D("h_invXsec_etaprime_etaLarge","#eta' in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_etaprime_etaTPC   = new TH1D("h_invXsec_etaprime_etaTPC"  ,"#eta' in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_etaprime_etaEMCal = new TH1D("h_invXsec_etaprime_etaEMCal","#eta' in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_etaprime_etaPHOS  = new TH1D("h_invXsec_etaprime_etaPHOS" ,"#eta' in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // omega meson
  TH1D *h_invXsec_omega_yDefault = new TH1D("h_invXsec_omega_yDefault","#omega in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_omega_etaLarge = new TH1D("h_invXsec_omega_etaLarge","#omega in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_omega_etaTPC   = new TH1D("h_invXsec_omega_etaTPC"  ,"#omega in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_omega_etaEMCal = new TH1D("h_invXsec_omega_etaEMCal","#omega in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_omega_etaPHOS  = new TH1D("h_invXsec_omega_etaPHOS" ,"#omega in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // direct photons (consider only direct photons)
  TH1D *h_invXsec_direct_photons_yDefault = new TH1D("h_invXsec_direct_photons_yDefault","direct photons in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_direct_photons_etaLarge = new TH1D("h_invXsec_direct_photons_etaLarge","direct photons in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_direct_photons_etaTPC   = new TH1D("h_invXsec_direct_photons_etaTPC"  ,"direct photons in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_direct_photons_etaEMCal = new TH1D("h_invXsec_direct_photons_etaEMCal","direct photons in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_direct_photons_etaPHOS  = new TH1D("h_invXsec_direct_photons_etaPHOS" ,"direct photons in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // shower/fragmentation photons only (gammas from "q -> q gamma" splitting)
  TH1D *h_invXsec_shower_photons_yDefault = new TH1D("h_invXsec_shower_photons_yDefault","shower photons (q -> q #gamma) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_shower_photons_etaLarge = new TH1D("h_invXsec_shower_photons_etaLarge","shower photons (q -> q #gamma) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_shower_photons_etaTPC   = new TH1D("h_invXsec_shower_photons_etaTPC"  ,"shower photons (q -> q #gamma) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_shower_photons_etaEMCal = new TH1D("h_invXsec_shower_photons_etaEMCal","shower photons (q -> q #gamma) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_shower_photons_etaPHOS  = new TH1D("h_invXsec_shower_photons_etaPHOS" ,"shower photons (q -> q #gamma) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // photons from ME (aka prompt)
  TH1D *h_invXsec_222_photons_yDefault = new TH1D("h_invXsec_222_photons_yDefault","photons from ME (aka prompt) in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_222_photons_etaLarge = new TH1D("h_invXsec_222_photons_etaLarge","photons from ME (aka prompt) in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_222_photons_etaTPC   = new TH1D("h_invXsec_222_photons_etaTPC"  ,"photons from ME (aka prompt) in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_222_photons_etaEMCal = new TH1D("h_invXsec_222_photons_etaEMCal","photons from ME (aka prompt) in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_222_photons_etaPHOS  = new TH1D("h_invXsec_222_photons_etaPHOS" ,"photons from ME (aka prompt) in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // decay photons
  TH1D *h_invXsec_decay_photons_yDefault = new TH1D("h_invXsec_decay_photons_yDefault","decay photons in |y| < 0.8", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_decay_photons_etaLarge = new TH1D("h_invXsec_decay_photons_etaLarge","decay photons in |#eta| < 3.0", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_decay_photons_etaTPC   = new TH1D("h_invXsec_decay_photons_etaTPC"  ,"decay photons in |#eta| < 0.9", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_decay_photons_etaEMCal = new TH1D("h_invXsec_decay_photons_etaEMCal","decay photons in |#eta| < 0.66", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_decay_photons_etaPHOS  = new TH1D("h_invXsec_decay_photons_etaPHOS" ,"decay photons in |#eta| < 0.12", pyHelp.ptBins, pyHelp.ptBinArray);

  // isolated photons (considers only direct photons)
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R03_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R03_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R04_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R04_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaTPC   = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaTPC","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaEMCal","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged2GeV_R05_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged2GeV_R05_photons_etaPHOS","direct iso (charged pt 2 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R03_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R03_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R04_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R04_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaTPC   = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaTPC","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaEMCal = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaEMCal","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_charged3GeV_R05_photons_etaPHOS  = new TH1D("h_invXsec_iso_charged3GeV_R05_photons_etaPHOS","direct iso (charged pt 3 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R03_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R03_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R04_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R04_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaTPC      = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaTPC","direct iso (full pt 2 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaEMCal    = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaEMCal","direct iso (full pt 2 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full2GeV_R05_photons_etaPHOS     = new TH1D("h_invXsec_iso_full2GeV_R05_photons_etaPHOS","direct iso (full pt 2 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.3) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.3) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R03_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R03_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.3) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.4) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.4) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R04_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R04_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.4) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);

  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaTPC      = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaTPC","direct iso (full pt 3 GeV/c in R=0.5) photons_etaTPC", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaEMCal    = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaEMCal","direct iso (full pt 3 GeV/c in R=0.5) photons_etaEMCal", pyHelp.ptBins, pyHelp.ptBinArray);
  TH1D *h_invXsec_iso_full3GeV_R05_photons_etaPHOS     = new TH1D("h_invXsec_iso_full3GeV_R05_photons_etaPHOS","direct iso (full pt 3 GeV/c in R=0.5) photons_etaPHOS", pyHelp.ptBins, pyHelp.ptBinArray);



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

  vector <TH1D*> vec_electron_yDefault_bin;
  vector <TH1D*> vec_electron_etaLarge_bin;
  vector <TH1D*> vec_electron_etaTPC_bin;
  vector <TH1D*> vec_electron_etaEMCal_bin;
  vector <TH1D*> vec_electron_etaPHOS_bin;
  
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

  //
  // ─── PLOTS FOR ISOLATION ITSELF ─────────────────────────────────────────────────
  //

  // decay photons

  vector<TH1D *> vec_iso_charged_R03_decay_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R04_decay_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R05_decay_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_charged_R03_decay_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R04_decay_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R05_decay_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_charged_R03_decay_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R04_decay_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R05_decay_photons_etaPHOS_bin;

  vector<TH1D *> vec_iso_full_R03_decay_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R04_decay_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R05_decay_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_full_R03_decay_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R04_decay_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R05_decay_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_full_R03_decay_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R04_decay_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R05_decay_photons_etaPHOS_bin;

  // direct photons
  vector<TH1D *> vec_iso_charged_R03_direct_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R04_direct_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R05_direct_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_charged_R03_direct_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R04_direct_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R05_direct_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_charged_R03_direct_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R04_direct_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R05_direct_photons_etaPHOS_bin;

  vector<TH1D *> vec_iso_full_R03_direct_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R04_direct_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R05_direct_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_full_R03_direct_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R04_direct_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R05_direct_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_full_R03_direct_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R04_direct_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R05_direct_photons_etaPHOS_bin;

  // all photons
  vector<TH1D *> vec_iso_charged_R03_all_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R04_all_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_charged_R05_all_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_charged_R03_all_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R04_all_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_charged_R05_all_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_charged_R03_all_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R04_all_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_charged_R05_all_photons_etaPHOS_bin;

  vector<TH1D *> vec_iso_full_R03_all_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R04_all_photons_etaTPC_bin;
  vector<TH1D *> vec_iso_full_R05_all_photons_etaTPC_bin;

  vector<TH1D *> vec_iso_full_R03_all_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R04_all_photons_etaEMCal_bin;
  vector<TH1D *> vec_iso_full_R05_all_photons_etaEMCal_bin;

  vector<TH1D *> vec_iso_full_R03_all_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R04_all_photons_etaPHOS_bin;
  vector<TH1D *> vec_iso_full_R05_all_photons_etaPHOS_bin;

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

    // electrons
    vec_electron_yDefault_bin.push_back( (TH1D*)h_electron_yDefault->Clone(Form("h_electron_yDefault_bin_%02d",i)) );
    vec_electron_etaLarge_bin.push_back( (TH1D*)h_electron_etaLarge->Clone(Form("h_electron_etaLarge_bin_%02d",i)) );
    vec_electron_etaTPC_bin.push_back( (TH1D*)h_electron_etaTPC->Clone(Form("h_electron_etaTPC_bin_%02d",i)) );
    vec_electron_etaEMCal_bin.push_back( (TH1D*)h_electron_etaEMCal->Clone(Form("h_electron_etaEMCal_bin_%02d",i)) );
    vec_electron_etaPHOS_bin.push_back( (TH1D*)h_electron_etaPHOS->Clone(Form("h_electron_etaPHOS_bin_%02d",i)) );

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

    if(producePhotonIsoSpectra){
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

      //
      // ─── SUM OF PT IN CONE ───────────────────────────────────────────
      //
      
      // decay photon histos TPC
      vec_iso_charged_R03_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R03_decay_photons_etaTPC->Clone(Form("h_iso_charged_R03_decay_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R04_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R04_decay_photons_etaTPC->Clone(Form("h_iso_charged_R04_decay_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R05_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R05_decay_photons_etaTPC->Clone(Form("h_iso_charged_R05_decay_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R03_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R03_decay_photons_etaTPC->Clone(Form("h_iso_full_R03_decay_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R04_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R04_decay_photons_etaTPC->Clone(Form("h_iso_full_R04_decay_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R05_decay_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R05_decay_photons_etaTPC->Clone(Form("h_iso_full_R05_decay_photons_etaTPC_bin_%02d", i)));

      // decay photon histos EMCAL
      vec_iso_charged_R03_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R03_decay_photons_etaEMCal->Clone(Form("h_iso_charged_R03_decay_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R04_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R04_decay_photons_etaEMCal->Clone(Form("h_iso_charged_R04_decay_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R05_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R05_decay_photons_etaEMCal->Clone(Form("h_iso_charged_R05_decay_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R03_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R03_decay_photons_etaEMCal->Clone(Form("h_iso_full_R03_decay_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R04_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R04_decay_photons_etaEMCal->Clone(Form("h_iso_full_R04_decay_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R05_decay_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R05_decay_photons_etaEMCal->Clone(Form("h_iso_full_R05_decay_photons_etaEMCal_bin_%02d", i)));

      // decay photon histos PHOS
      vec_iso_charged_R03_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R03_decay_photons_etaPHOS->Clone(Form("h_iso_charged_R03_decay_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R04_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R04_decay_photons_etaPHOS->Clone(Form("h_iso_charged_R04_decay_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R05_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R05_decay_photons_etaPHOS->Clone(Form("h_iso_charged_R05_decay_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R03_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R03_decay_photons_etaPHOS->Clone(Form("h_iso_full_R03_decay_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R04_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R04_decay_photons_etaPHOS->Clone(Form("h_iso_full_R04_decay_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R05_decay_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R05_decay_photons_etaPHOS->Clone(Form("h_iso_full_R05_decay_photons_etaPHOS_bin_%02d", i)));

      // direct photon histos TPC
      vec_iso_charged_R03_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R03_direct_photons_etaTPC->Clone(Form("h_iso_charged_R03_direct_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R04_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R04_direct_photons_etaTPC->Clone(Form("h_iso_charged_R04_direct_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R05_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R05_direct_photons_etaTPC->Clone(Form("h_iso_charged_R05_direct_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R03_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R03_direct_photons_etaTPC->Clone(Form("h_iso_full_R03_direct_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R04_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R04_direct_photons_etaTPC->Clone(Form("h_iso_full_R04_direct_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R05_direct_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R05_direct_photons_etaTPC->Clone(Form("h_iso_full_R05_direct_photons_etaTPC_bin_%02d", i)));

      // direct photon histos EMCAL
      vec_iso_charged_R03_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R03_direct_photons_etaEMCal->Clone(Form("h_iso_charged_R03_direct_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R04_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R04_direct_photons_etaEMCal->Clone(Form("h_iso_charged_R04_direct_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R05_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R05_direct_photons_etaEMCal->Clone(Form("h_iso_charged_R05_direct_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R03_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R03_direct_photons_etaEMCal->Clone(Form("h_iso_full_R03_direct_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R04_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R04_direct_photons_etaEMCal->Clone(Form("h_iso_full_R04_direct_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R05_direct_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R05_direct_photons_etaEMCal->Clone(Form("h_iso_full_R05_direct_photons_etaEMCal_bin_%02d", i)));

      // direct photon histos PHOS
      vec_iso_charged_R03_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R03_direct_photons_etaPHOS->Clone(Form("h_iso_charged_R03_direct_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R04_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R04_direct_photons_etaPHOS->Clone(Form("h_iso_charged_R04_direct_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R05_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R05_direct_photons_etaPHOS->Clone(Form("h_iso_charged_R05_direct_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R03_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R03_direct_photons_etaPHOS->Clone(Form("h_iso_full_R03_direct_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R04_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R04_direct_photons_etaPHOS->Clone(Form("h_iso_full_R04_direct_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R05_direct_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R05_direct_photons_etaPHOS->Clone(Form("h_iso_full_R05_direct_photons_etaPHOS_bin_%02d", i)));

      // all photon histos TPC
      vec_iso_charged_R03_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R03_all_photons_etaTPC->Clone(Form("h_iso_charged_R03_all_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R04_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R04_all_photons_etaTPC->Clone(Form("h_iso_charged_R04_all_photons_etaTPC_bin_%02d", i)));
      vec_iso_charged_R05_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_charged_R05_all_photons_etaTPC->Clone(Form("h_iso_charged_R05_all_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R03_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R03_all_photons_etaTPC->Clone(Form("h_iso_full_R03_all_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R04_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R04_all_photons_etaTPC->Clone(Form("h_iso_full_R04_all_photons_etaTPC_bin_%02d", i)));
      vec_iso_full_R05_all_photons_etaTPC_bin.push_back((TH1D *)h_iso_full_R05_all_photons_etaTPC->Clone(Form("h_iso_full_R05_all_photons_etaTPC_bin_%02d", i)));

      // all photon histos EMCAL
      vec_iso_charged_R03_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R03_all_photons_etaEMCal->Clone(Form("h_iso_charged_R03_all_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R04_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R04_all_photons_etaEMCal->Clone(Form("h_iso_charged_R04_all_photons_etaEMCal_bin_%02d", i)));
      vec_iso_charged_R05_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_charged_R05_all_photons_etaEMCal->Clone(Form("h_iso_charged_R05_all_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R03_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R03_all_photons_etaEMCal->Clone(Form("h_iso_full_R03_all_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R04_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R04_all_photons_etaEMCal->Clone(Form("h_iso_full_R04_all_photons_etaEMCal_bin_%02d", i)));
      vec_iso_full_R05_all_photons_etaEMCal_bin.push_back((TH1D *)h_iso_full_R05_all_photons_etaEMCal->Clone(Form("h_iso_full_R05_all_photons_etaEMCal_bin_%02d", i)));

      // all photon histos PHOS
      vec_iso_charged_R03_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R03_all_photons_etaPHOS->Clone(Form("h_iso_charged_R03_all_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R04_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R04_all_photons_etaPHOS->Clone(Form("h_iso_charged_R04_all_photons_etaPHOS_bin_%02d", i)));
      vec_iso_charged_R05_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_charged_R05_all_photons_etaPHOS->Clone(Form("h_iso_charged_R05_all_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R03_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R03_all_photons_etaPHOS->Clone(Form("h_iso_full_R03_all_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R04_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R04_all_photons_etaPHOS->Clone(Form("h_iso_full_R04_all_photons_etaPHOS_bin_%02d", i)));
      vec_iso_full_R05_all_photons_etaPHOS_bin.push_back((TH1D *)h_iso_full_R05_all_photons_etaPHOS->Clone(Form("h_iso_full_R05_all_photons_etaPHOS_bin_%02d", i)));
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


    if(producePhotonIsoSpectra){
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
      if(useGammaJetCorrelations){
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
	      TMath::Abs(p.event[i].eta()) < (etaTPC-jetRadius)){       // in maximal TPC-minus-iso-cone-radius acceptance

	    // photon as pseudojet for analysis
	    PseudoJet photonJet(p.event[i].px(), p.event[i].py(), p.event[i].pz(), p.event[i].e());
	    if(photonJet.pt() < 15.) continue;
	    if(photonJet.pt() > 30.) continue;
	    // calculate ue pt density for a given photon i NOT IMPLEMENTED IN THE MOMENT
	    double UEPtDensity = 0.;
	    // printf("UEPtDensity(p.event, i) = %f\n",UEPtDensity);
	    // check isolation
	    isPhotonIsolated = pyHelp.IsPhotonIsolated(p.event, i, etaTPC-jetRadius, isoConeRadius, isoPtMax, UEPtDensity, vec_isoCone_track_dPhi_bin.at(iBin), vec_isoCone_track_dEta_bin.at(iBin), vec_isoPt_bin.at(iBin), vec_isoPt_corrected_bin.at(iBin));

	    if(vJets.size() > 0 && isPhotonIsolated)
	      for(unsigned int iJet = 0; iJet < vJets.size(); iJet++){
		bool isJetSeparated = ( TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))) > TMath::Pi()/2. );
		if(vJets.at(iJet).pt() < 10.) break; // vJets are sorted by pt, break is ok
		if(iJet == 0) vec_xSecTriggerGamma_bin.at(iBin)->Fill(0.);  // if there is at least one jet, count trigger photons, but only once for all jets connected to this photon; can be used to normalize histograms per trigger photon in the end
		// gamma-jet correlation	 
		vec_dPhiJetGamma_bin.at(iBin)->Fill(TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))));
		if(!isJetSeparated) continue;
		// x_Jet-gamma
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
	delete cs;           
      } // end of "if(useGammaJetCorr)
      //------------------------------------------------------------------------------------------
      //----- END OF jets + photon correlation ---------------------------------------------------
      //------------------------------------------------------------------------------------------

      //------------------------------------------------------------------------------------------
      pyHelp.Fill_TH2_Electron_TopMotherID(p.event, etaEMCal, vec_electron_pt_topMotherID_bin.at(iBin));

      pyHelp.Fill_Electron_Pt(p.event, yDefault, true, vec_electron_yDefault_bin.at(iBin));
      pyHelp.Fill_Electron_Pt(p.event, etaLarge, false, vec_electron_etaLarge_bin.at(iBin));
      pyHelp.Fill_Electron_Pt(p.event, etaTPC, false, vec_electron_etaTPC_bin.at(iBin));
      pyHelp.Fill_Electron_Pt(p.event, etaEMCal, false, vec_electron_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Electron_Pt(p.event, etaPHOS, false, vec_electron_etaPHOS_bin.at(iBin));
      
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

      if(producePhotonIsoSpectra){
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


  //
  // ─── FILL SUM IN CONE ───────────────────────────────────────────────────────────
  //

  // fill isolated photons: considers only direct photons
  // arguments = (p.event, etaAcc, vec_histo, bool onlyCharged?, iso cone radius, iso pt)
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_charged_R03_decay_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_charged_R04_decay_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_charged_R05_decay_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R03_decay_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R04_decay_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R05_decay_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R03_decay_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R04_decay_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R05_decay_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);


  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_full_R03_decay_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_full_R04_decay_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaTPC, vec_iso_full_R05_decay_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_full_R03_decay_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_full_R04_decay_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaEMCal, vec_iso_full_R05_decay_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_full_R03_decay_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_full_R04_decay_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Decay_Photon_Pt(p.event, etaPHOS, vec_iso_full_R05_decay_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);

  //direct

  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_charged_R03_direct_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_charged_R04_direct_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_charged_R05_direct_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R03_direct_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R04_direct_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R05_direct_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R03_direct_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R04_direct_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R05_direct_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);

  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_full_R03_direct_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_full_R04_direct_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaTPC, vec_iso_full_R05_direct_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_full_R03_direct_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_full_R04_direct_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaEMCal, vec_iso_full_R05_direct_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_full_R03_direct_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_full_R04_direct_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_Direct_Photon_Pt(p.event, etaPHOS, vec_iso_full_R05_direct_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);

 // all
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_charged_R03_all_photons_etaTPC_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_charged_R04_all_photons_etaTPC_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_charged_R05_all_photons_etaTPC_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R03_all_photons_etaEMCal_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R04_all_photons_etaEMCal_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_charged_R05_all_photons_etaEMCal_bin.at(iBin), true, 0.5, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R03_all_photons_etaPHOS_bin.at(iBin), true, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R04_all_photons_etaPHOS_bin.at(iBin), true, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_charged_R05_all_photons_etaPHOS_bin.at(iBin), true, 0.5, 2.);

  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_full_R03_all_photons_etaTPC_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_full_R04_all_photons_etaTPC_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaTPC, vec_iso_full_R05_all_photons_etaTPC_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_full_R03_all_photons_etaEMCal_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_full_R04_all_photons_etaEMCal_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaEMCal, vec_iso_full_R05_all_photons_etaEMCal_bin.at(iBin), false, 0.5, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_full_R03_all_photons_etaPHOS_bin.at(iBin), false, 0.3, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_full_R04_all_photons_etaPHOS_bin.at(iBin), false, 0.4, 2.);
  pyHelp.Fill_iso_All_Photon_Pt(p.event, etaPHOS, vec_iso_full_R05_all_photons_etaPHOS_bin.at(iBin), false, 0.5, 2.);
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

      if(producePhotonIsoSpectra){
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
    cout << "sigma = " << sigma << endl;
    cout << "weightSum = " << p.info.weightSum() << endl;

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

    vec_electron_yDefault_bin.at(iBin)->Scale(sigma);
    vec_electron_etaLarge_bin.at(iBin)->Scale(sigma);
    vec_electron_etaTPC_bin.at(iBin)->Scale(sigma);
    vec_electron_etaEMCal_bin.at(iBin)->Scale(sigma);
    vec_electron_etaPHOS_bin.at(iBin)->Scale(sigma);
  
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

    if(producePhotonIsoSpectra){
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

      //
      // ─── PT IN CONE ──────────────────────────────────────────────────
      //

      vec_iso_charged_R03_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);


      vec_iso_full_R03_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      // direct
      vec_iso_charged_R03_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      vec_iso_full_R03_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_direct_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_direct_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_direct_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      // all
      vec_iso_charged_R03_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R03_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R04_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_charged_R05_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);

      vec_iso_full_R03_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_all_photons_etaTPC_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_all_photons_etaEMCal_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R03_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R04_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
      vec_iso_full_R05_all_photons_etaPHOS_bin.at(iBin)->Scale(sigma);
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


    if(producePhotonIsoSpectra){
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
  // write info string to file
  pyHelp.Write_README(p, file, argc, argv, pdfA, pdfB, pTHatBin);

  // store weightSum for normalization later on (in standalone Pythia = number of events)
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
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_yDefault_bin, h_electron_yDefault, file, dir_electron, 2*yDefault, true);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_etaLarge_bin, h_electron_etaLarge, file, dir_electron, 2*etaLarge, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_etaTPC_bin, h_electron_etaTPC, file, dir_electron, 2*etaTPC, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_etaEMCal_bin, h_electron_etaEMCal, file, dir_electron, 2*etaEMCal, false);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_etaPHOS_bin, h_electron_etaPHOS, file, dir_electron, 2*etaPHOS, false);

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

  if(producePhotonIsoSpectra){
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


    //
    // ─── PT IN CONE ──────────────────────────────────────────────────
    //

    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_decay_photons_etaTPC_bin, h_iso_charged_R03_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_decay_photons_etaTPC_bin, h_iso_charged_R04_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_decay_photons_etaTPC_bin, h_iso_charged_R05_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_decay_photons_etaEMCal_bin, h_iso_charged_R03_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_decay_photons_etaEMCal_bin, h_iso_charged_R04_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_decay_photons_etaEMCal_bin, h_iso_charged_R05_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_decay_photons_etaPHOS_bin, h_iso_charged_R03_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_decay_photons_etaPHOS_bin, h_iso_charged_R04_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_decay_photons_etaPHOS_bin, h_iso_charged_R05_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);

    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_decay_photons_etaTPC_bin, h_iso_full_R03_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_decay_photons_etaTPC_bin, h_iso_full_R04_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_decay_photons_etaTPC_bin, h_iso_full_R05_decay_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_decay_photons_etaEMCal_bin, h_iso_full_R03_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_decay_photons_etaEMCal_bin, h_iso_full_R04_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_decay_photons_etaEMCal_bin, h_iso_full_R05_decay_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_decay_photons_etaPHOS_bin, h_iso_full_R03_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_decay_photons_etaPHOS_bin, h_iso_full_R04_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_decay_photons_etaPHOS_bin, h_iso_full_R05_decay_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);

    // direct
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_direct_photons_etaTPC_bin, h_iso_charged_R03_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_direct_photons_etaTPC_bin, h_iso_charged_R04_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_direct_photons_etaTPC_bin, h_iso_charged_R05_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_direct_photons_etaEMCal_bin, h_iso_charged_R03_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_direct_photons_etaEMCal_bin, h_iso_charged_R04_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_direct_photons_etaEMCal_bin, h_iso_charged_R05_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_direct_photons_etaPHOS_bin, h_iso_charged_R03_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_direct_photons_etaPHOS_bin, h_iso_charged_R04_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_direct_photons_etaPHOS_bin, h_iso_charged_R05_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);

    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_direct_photons_etaTPC_bin, h_iso_full_R03_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_direct_photons_etaTPC_bin, h_iso_full_R04_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_direct_photons_etaTPC_bin, h_iso_full_R05_direct_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_direct_photons_etaEMCal_bin, h_iso_full_R03_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_direct_photons_etaEMCal_bin, h_iso_full_R04_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_direct_photons_etaEMCal_bin, h_iso_full_R05_direct_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_direct_photons_etaPHOS_bin, h_iso_full_R03_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_direct_photons_etaPHOS_bin, h_iso_full_R04_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_direct_photons_etaPHOS_bin, h_iso_full_R05_direct_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    
    // all
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_all_photons_etaTPC_bin, h_iso_charged_R03_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_all_photons_etaTPC_bin, h_iso_charged_R04_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_all_photons_etaTPC_bin, h_iso_charged_R05_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_all_photons_etaEMCal_bin, h_iso_charged_R03_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_all_photons_etaEMCal_bin, h_iso_charged_R04_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_all_photons_etaEMCal_bin, h_iso_charged_R05_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R03_all_photons_etaPHOS_bin, h_iso_charged_R03_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R04_all_photons_etaPHOS_bin, h_iso_charged_R04_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_charged_R05_all_photons_etaPHOS_bin, h_iso_charged_R05_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);

    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_all_photons_etaTPC_bin, h_iso_full_R03_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_all_photons_etaTPC_bin, h_iso_full_R04_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_all_photons_etaTPC_bin, h_iso_full_R05_all_photons_etaTPC, file, dir_isoGamma, 2 * etaTPC, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_all_photons_etaEMCal_bin, h_iso_full_R03_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_all_photons_etaEMCal_bin, h_iso_full_R04_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_all_photons_etaEMCal_bin, h_iso_full_R05_all_photons_etaEMCal, file, dir_isoGamma, 2 * etaEMCal, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R03_all_photons_etaPHOS_bin, h_iso_full_R03_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R04_all_photons_etaPHOS_bin, h_iso_full_R04_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
    pyHelp.Add_Histos_Scale_Write2File(vec_iso_full_R05_all_photons_etaPHOS_bin, h_iso_full_R05_all_photons_etaPHOS, file, dir_isoGamma, 2 * etaPHOS, false);
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

  if(producePhotonIsoSpectra){
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
