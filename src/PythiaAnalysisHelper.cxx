#include "PythiaAnalysisHelper.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "fastjet/ClusterSequence.hh"

using fastjet::PseudoJet;
using fastjet::JetDefinition;
using fastjet::ClusterSequence;
using fastjet::antikt_algorithm;


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Set_Pythia_Randomseed(Pythia8::Pythia &p){

  TRandom3 *rand = new TRandom3(0);
  p.readString("Random:setSeed = on");
  p.readString(Form("Random:seed = %ld", (long)(rand->Rndm()*9e8) ));
  delete rand;
  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv){

  p.readString(Form("Beams:eCM = %f", strtof(argv[4], NULL) ));
  if(argc > 6)
    p.readString(Form("SigmaProcess:renormMultFac = %f", strtof(argv[6],NULL) ));
  if(argc > 7){
    p.readString(Form("SigmaProcess:factorMultFac = %f", strtof(argv[7],NULL) ));
  }

  if(argc > 5){
    if( !strcmp(argv[5],"noHadro") ){
      p.readString("HadronLevel:all = off");
      printf("\nPythia events are generated without Hadronization\n");
    }
    else if( !strcmp(argv[5],"noMPI") ){
      p.readString("PartonLevel:MPI = off");
      printf("\nPythia events are generated without Multiparton Interaction\n");
    }
    else if( !strcmp(argv[5],"noMPInoHadro") ){
      p.readString("PartonLevel:MPI = off");
      p.readString("HadronLevel:all = off");
      printf("\nPythia events are generated without Multiparton Interaction and Hadronization\n");
    }
    else if( !strcmp(argv[5],"noShower") ){
      p.readString("Check:event = off");
      p.readString("PartonLevel:MPI = off");
      p.readString("PartonLevel:FSR = off");
      p.readString("PartonLevel:ISR = off");
      p.readString("PartonLevel:Remnants = off");
      p.readString("HadronLevel:all = off");
      printf("\nPythia events are generated without parton shower (only LO scattering)\n");
    }
    else if( !strcmp(argv[5],"fullEvents") ){
      printf("\nFull pythia events will be generated (default)\n");
    }
    else if (!strcmp(argv[5], "fullEventsMonash"))
    {
      p.readString("Tune:pp = 14");
      printf("\nFull pythia events will be generated using Monash 2013 tune\n");
    }
    else if (!strcmp(argv[5], "fullEventsTune4C"))
    {
      p.readString("Tune:pp = 5");
      printf("\nFull pythia events will be generated using Tune 4C tune\n");
    }
    else
      printf("\nNo sensible argument argv[5] is given -> full pythia events will be generated (default)\n");
  }
  
  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Write_README(Pythia8::Pythia &p, TFile &file, int argc, char **argv, string pdfA, string pdfB, double *pTHatBin){

  file.cd();
  string infoPDFA, infoPDFB, infoEnergyA, infoEnergyB, infoProcess, infoOption, infoRenScale, infoFacScale;

  if( !strcmp(argv[2],"JJ") ) infoProcess = "Using JetJet processes (codes 111-124)";
  else if( !strcmp(argv[2],"GJ") ) infoProcess = "Using GammaJet processes (codes 201-205)";
  else if( !strcmp(argv[2],"WeakBoson" ) ) infoProcess = "Using single weak boson production processes (codes 224,225)";
  else if( !strcmp(argv[2],"MB" ) ) infoProcess = "Using Minimum Bias production processes (codes 101,103-106)";
  else if( !strcmp(argv[2],"MBVeto" ) ) infoProcess = Form("Using Minimum Bias production processes (codes 101,103-106), but only below pthat = %f GeV, so that other processes can be used for high pt statistics",pTHatBin[0]);

  infoEnergyA = Form("%f", p.event[1].e());
  infoEnergyB = Form("%f", p.event[2].e());

  if( !strcmp(argv[5],"fullEvents") ) infoOption = "Full events have been processed (including hadronization, MPI etc...)";
  if( !strcmp(argv[5],"fullEventsMonash") ) infoOption = "Full events have been processed (including hadronization, MPI etc...) using the Monash 2013 tune";
  if( !strcmp(argv[5],"fullEventsTune4C") ) infoOption = "Full events have been processed (including hadronization, MPI etc...) using the Tune 4C";
  if( !strcmp(argv[5],"noMPI") ) infoOption = "BEWARE: events have been processed without Multiparton Interaction";
  if( !strcmp(argv[5],"noMPInoHadro") ) infoOption = "BEWARE: events have been processed without Multiparton Interaction and without hadronization";
  if( !strcmp(argv[5],"noShower") ) infoOption = "BEWARE: events have been processed without any parton shower, only the naked hard process";
    
  if( argc >= 7 ) infoRenScale = Form("Using renormalization scale: #mu_{R} = Q^{2} #times %s", argv[6]);
  else infoRenScale = "Using renormalization scale: #mu^{2}_{R} = Q^{2} #times 1.00";
  if( argc >= 8 ) infoFacScale = Form("Using factorization scale: #mu_{F} = Q^{2} #times %s", argv[7]);
  else infoFacScale = "Using factorization scale: #mu^{2}_{F} = Q^{2} #times 1.00";

  if(argc >= 10) infoPDFA = pdfA;
  else infoPDFA = "NNPDF2.3 QCD+QED LO";

  if(argc >= 11) infoPDFB = pdfB; 
  else infoPDFB = "NNPDF2.3 QCD+QED LO";

  TCanvas *c = new TCanvas(Form("README_%s",argv[2]),Form("README_%s",argv[2]));
  TPaveText *pt = new TPaveText(.01,.01,.99,.99);
  pt->AddText("Generated with Pythia 8.230 and following config:");
  pt->AddText(Form("PDF beam A = %s", infoPDFA.c_str()));
  pt->AddText(Form("PDF beam B = %s", infoPDFB.c_str()));
  pt->AddText(Form("Energy beam A = %s GeV", infoEnergyA.c_str()));
  pt->AddText(Form("Energy beam B = %s GeV", infoEnergyB.c_str()));
  pt->AddText(Form("%s", infoOption.c_str()));
  pt->AddText(Form("%s", infoRenScale.c_str()));
  pt->AddText(Form("%s", infoFacScale.c_str()));
  pt->AddText(Form("%s", infoProcess.c_str()));
  pt->AddText("Check other READMEs for additional processes that may have been included by merging.");
  pt->Draw();
  c->Write();

  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::ProcessSwitch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p){

  if( !strcmp(argv[2],"MB") || !strcmp(argv[2],"MBVeto")){
    //p.readString("SoftQCD:inelastic= on"); // includes also double- and single diffractive events, only necessary for lowest pt
    p.readString("SoftQCD:nonDiffractive= on");
  return;
  }
  else if( !strcmp(argv[2],"JJ") ){
    p.readString("HardQCD:all = on");
    p.readString("HardQCD:hardccbar = on");
    p.readString("HardQCD:hardbbbar = on");
  }
  else if( !strcmp(argv[2],"GJ") ){
    p.readString("PromptPhoton:all = on");
  }
  else if( !strcmp(argv[2],"WeakBoson") ){
    // the two commented processes cannot be used with pthat bins
    // p.readString("WeakBosonExchange:all = on");
    // p.readString("WeakSingleBoson:all = on");
    // instead, use the following
    // note that there are more weak boson processes, e.g. double boson production
    // or fermion production via t-channel weak boson exchange
    // but the cross section is very small compared to single boson production
    p.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on");
    p.readString("WeakSingleBoson:ffbar2ffbar(s:W) = on");
  }
  else
    printf("\nNo process switched on. Provide \"MB\" or \"MBVeto\" or \"JJ\" or \"GJ\" or \"WeakBoson\" as second argument\n");
  
  p.settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
  p.settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);

  return;
}

//----------------------------------------------------------------------
// PYTHIA8's phi goes from -pi to pi; compute correct angle difference
double PythiaAnalysisHelper::CorrectPhiDelta(double angle1, double angle2){
  double pi = TMath::Pi();
  double phi = TMath::Abs(angle1 - angle2);
  if(phi >= pi) return 2*pi-phi;
  else return phi;
}

//----------------------------------------------------------------------
// p going in +z direction
double PythiaAnalysisHelper::XObs_pGoing(PseudoJet &hadronjet, PseudoJet &photonjet, double beamEnergy){
  double xobs;
  xobs = (hadronjet.pt()*TMath::Exp(-hadronjet.eta()) + photonjet.pt()*TMath::Exp(-photonjet.eta())) /(2*beamEnergy);
  return xobs;
}

//----------------------------------------------------------------------
// Pb going in +z direction
double PythiaAnalysisHelper::XObs_PbGoing(PseudoJet &hadronjet, PseudoJet &photonjet, double beamEnergy){
  double xobs;
  //xobs = (hadronjet.pt()*TMath::Exp(-hadronjet.eta()) + photonjet.pt()*TMath::Exp(-photonjet.eta())) /(2*beamEnergy);
  xobs = (hadronjet.pt()*TMath::Exp(hadronjet.eta()) + photonjet.pt()*TMath::Exp(photonjet.eta())) /(2*beamEnergy);
  return xobs;
}

//----------------------------------------------------------------------
// isolation cut: sum energy around photon and abandon event if threshold is reached
bool PythiaAnalysisHelper::IsPhotonIsolated(Event &event, int iPhoton, const double &etaAbsMaxPhoton, const double &isoConeRadius, const double &isoPtMax, double UEPtDensity, TH1D *h_phi, TH1D *h_eta, TH1D *h_isoPt, TH1D *h_isoPt_corrected){

  double isoCone_dR = 999.;
  double isoCone_pt = 0.; // reset sum of energy in cone
  
  for (int iTrack = 5; iTrack < event.size(); iTrack++) {
    if ( !event[iTrack].isFinal() ) continue;
    if ( !event[iTrack].isVisible() ) continue;
    if ( !event[iTrack].isCharged() ) continue;
    if ( TMath::Abs(event[iPhoton].eta()) > etaAbsMaxPhoton ) continue;
    if ( iTrack == iPhoton ) continue; // dont count photon, not necessary for tracks obviously

    // distance between photon and particle at index iTrack
    isoCone_dR = sqrt( pow(CorrectPhiDelta(event[iTrack].phi(), event[iPhoton].phi()), 2)
		       + pow(event[iTrack].eta() - event[iPhoton].eta(), 2) );
	    
    if(isoCone_dR < isoConeRadius){
      isoCone_pt += event[iTrack].pT();
      h_phi->Fill(event[iTrack].phi() - event[iPhoton].phi());
      h_eta->Fill(event[iTrack].eta() - event[iPhoton].eta());
    }
    h_isoPt_corrected->Fill(isoCone_pt-(UEPtDensity*0.4*0.4*TMath::Pi()));
    h_isoPt->Fill(isoCone_pt);
  }
      
  if( isoCone_pt >= isoPtMax ) return false;
  else return true;
}
//----------------------------------------------------------------------
bool PythiaAnalysisHelper::IsPhotonIsolatedPowheg(Event &event, int iPhoton, const double &etaAbsMaxPhoton, const double &isoConeRadius, const double &isoPtMax, double UEPtDensity, vector<TH1D> &vec_phi, vector<TH1D> &vec_eta, vector<TH1D> &vec_isoPt, vector<TH1D> &vec_isoPt_corrected, vector<double> vec_weights){

  double isoCone_dR = 999.;
  double isoCone_pt = 0.; // reset sum of energy in cone
  
  for (int iTrack = 5; iTrack < event.size(); iTrack++) {
    if ( !event[iTrack].isFinal() ) continue;
    if ( !event[iTrack].isVisible() ) continue;
    if ( !event[iTrack].isCharged() ) continue;
    if ( TMath::Abs(event[iPhoton].eta()) > etaAbsMaxPhoton ) continue;
    if ( iTrack == iPhoton ) continue; // dont count photon, not necessary for tracks obviously

    // distance between photon and particle at index iTrack
    isoCone_dR = sqrt( pow(CorrectPhiDelta(event[iTrack].phi(), event[iPhoton].phi()), 2)
		       + pow(event[iTrack].eta() - event[iPhoton].eta(), 2) );
	    
    if(isoCone_dR < isoConeRadius){
      isoCone_pt += event[iTrack].pT();
      FillForEachWeight(vec_phi, event[iTrack].phi() - event[iPhoton].phi(), vec_weights);
      FillForEachWeight(vec_eta, event[iTrack].eta() - event[iPhoton].eta(), vec_weights);
    }
    FillForEachWeight(vec_isoPt_corrected, isoCone_pt-(UEPtDensity*0.4*0.4*TMath::Pi()), vec_weights);
    FillForEachWeight(vec_isoPt, isoCone_pt, vec_weights);
  }
      
  if( isoCone_pt >= isoPtMax ) return false;
  else return true;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Electron_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT());
    }else{
      if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT());      
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Pi0_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 111 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT());
    }else{
      if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT());      
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 111 && TMath::Abs(event[i].y()) < etaMax ) {
	int mI = event[i].mother1();
	if ( !(TMath::Abs(event[mI].id()) == 310   || // K0_s, K0_l
	       TMath::Abs(event[mI].id()) == 321   || // K+,K-
	       TMath::Abs(event[mI].id()) == 3122  || // Lambda, Anti-Lambda
	       TMath::Abs(event[mI].id()) == 3212  || // Sigma0
	       TMath::Abs(event[mI].id()) == 3222  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3112  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3322  || // Cascades
	       TMath::Abs(event[mI].id()) == 3312)  ) // Cascades
	  {
	    h->Fill(event[i].pT());
	  }
      }
    }else{
      if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) {
	int mI = event[i].mother1();
	if ( !(TMath::Abs(event[mI].id()) == 310   || // K0_s, K0_l
	       TMath::Abs(event[mI].id()) == 321   || // K+,K-
	       TMath::Abs(event[mI].id()) == 3122  || // Lambda, Anti-Lambda
	       TMath::Abs(event[mI].id()) == 3212  || // Sigma0
	       TMath::Abs(event[mI].id()) == 3222  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3112  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3322  || // Cascades
	       TMath::Abs(event[mI].id()) == 3312)  ) // Cascades
	  {
	    h->Fill(event[i].pT());
	  }
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Eta_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 221 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT());
    }else{
      if(event[i].id() == 221 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT());      
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_EtaPrime_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 331 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT());
    }else{
      if(event[i].id() == 331 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT());      
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Omega_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 223 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT());
    }else{
      if(event[i].id() == 223 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT());      
    } 
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].y()) < etaMax && event[i].status() < 90 )
      h->Fill(event[i].pT()); 
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax && event[i].status() < 90){
	if( TMath::Abs(event[event[i].iTopCopy()].status() ) > 40 )
	  h->Fill(event[i].pT());
      }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90)
	if( TMath::Abs(event[event[i].iTopCopy()].status() ) < 40 )
	  h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Electron_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_ElectronNeg_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_ElectronPos_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == -11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                                  bool isoCharged, double iso_cone_radius, double iso_pt){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90){
        // isolation check------------------------------
        double pt_temp = 0.;
        if( isoCharged )
          for(int j = 5; j < event.size(); j++){
            // only charged considered for iso cut
            if( event[j].isFinal() && event[j].isVisible() && event[j].isCharged() && j != i)
              if( TMath::Sqrt(   (event[i].phi()-event[j].phi()) * (event[i].phi()-event[j].phi())
                                 + (event[i].eta()-event[j].eta()) * (event[i].eta()-event[j].eta()) )
                  < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        // charged + neutral pt considered for iso cut
        if( !isoCharged )
          for(int j = 5; j < event.size(); j++){
            if( event[j].isFinal() && event[j].isVisible() && j != i)
              if( TMath::Sqrt(   (event[i].phi()-event[j].phi()) * (event[i].phi()-event[j].phi())
                                 + (event[i].eta()-event[j].eta()) * (event[i].eta()-event[j].eta()) )
                  < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        if( pt_temp <= iso_pt)
          h->Fill(event[i].pT());

        //----------------------------------------------
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() > 90) h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Pi0_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 111 && TMath::Abs(event[i].y()) < etaMax ) h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }else{
      if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 111 && TMath::Abs(event[i].y()) < etaMax ) {
	int mI = event[i].mother1();
	if ( !(TMath::Abs(event[mI].id()) == 310   || // K0_s, K0_l
	       TMath::Abs(event[mI].id()) == 321   || // K+,K-
	       TMath::Abs(event[mI].id()) == 3122  || // Lambda, Anti-Lambda
	       TMath::Abs(event[mI].id()) == 3212  || // Sigma0
	       TMath::Abs(event[mI].id()) == 3222  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3112  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3322  || // Cascades
	       TMath::Abs(event[mI].id()) == 3312)  ) // Cascades
	  {
	    h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
	  }
      }
    }else{
      if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) {
	int mI = event[i].mother1();
	if ( !(TMath::Abs(event[mI].id()) == 310   || // K0_s, K0_l
	       TMath::Abs(event[mI].id()) == 321   || // K+,K-
	       TMath::Abs(event[mI].id()) == 3122  || // Lambda, Anti-Lambda
	       TMath::Abs(event[mI].id()) == 3212  || // Sigma0
	       TMath::Abs(event[mI].id()) == 3222  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3112  || // Sigmas
	       TMath::Abs(event[mI].id()) == 3322  || // Cascades
	       TMath::Abs(event[mI].id()) == 3312)  ) // Cascades
	  {
	    h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
	  }
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Eta_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 221 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
    }else{
      if(event[i].id() == 221 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));      
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_EtaPrime_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 331 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
    }else{
      if(event[i].id() == 331 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Omega_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(useRap){
      if(event[i].id() == 223 && TMath::Abs(event[i].y()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));
    }else{
      if(event[i].id() == 223 && TMath::Abs(event[i].eta()) < etaMax ) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()));      
    } 
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90)
	if( TMath::Abs(event[event[i].iTopCopy()].status() ) > 40 )
	  h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90)
	if( TMath::Abs(event[event[i].iTopCopy()].status() ) < 40 )
	  h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                                           bool isoCharged, double iso_cone_radius, double iso_pt){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90){
        // isolation check------------------------------
        double pt_temp = 0.;
        if( isoCharged )
          for(int j = 5; j < event.size(); j++){
            // only charged considered for iso cut
            if( event[j].isFinal() && event[j].isVisible() && event[j].isCharged() && j != i)
              if( TMath::Sqrt(   (event[i].phi()-event[j].phi()) * (event[i].phi()-event[j].phi())
                                 + (event[i].eta()-event[j].eta()) * (event[i].eta()-event[j].eta()) )
                  < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        // charged + neutral pt considered for iso cut
        if( !isoCharged )
          for(int j = 5; j < event.size(); j++){
            if( event[j].isFinal() && event[j].isVisible() && j != i)
              if( TMath::Sqrt(   (event[i].phi()-event[j].phi()) * (event[i].phi()-event[j].phi())
                                 + (event[i].eta()-event[j].eta()) * (event[i].eta()-event[j].eta()) )
                  < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        if( pt_temp <= iso_pt)
          h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );

        //----------------------------------------------
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_iso_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                                             bool isoCharged, double iso_cone_radius, double iso_pt)
{
  for (int i = 5; i < event.size(); i++)
  {
    if (event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax)
    {
      if (event[i].status() < 90)
      {
        // isolation check------------------------------
        double pt_temp = 0.;
        if (isoCharged)
          for (int j = 5; j < event.size(); j++)
          {
            // only charged considered for iso cut
            if (event[j].isFinal() && event[j].isVisible() && event[j].isCharged() && j != i)
              if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        // charged + neutral pt considered for iso cut
        if (!isoCharged)
          for (int j = 5; j < event.size(); j++)
          {
            if (event[j].isFinal() && event[j].isVisible() && j != i)
              if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
                pt_temp += event[j].pT();
          }

          h->Fill(pt_temp);
        //----------------------------------------------
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_iso_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                                         bool isoCharged, double iso_cone_radius, double iso_pt)
{
  for (int i = 5; i < event.size(); i++)
  {
    if (event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax)
    {
      if ((event[i].status() > 90) && (event[i].status() < 100))
      {
        // isolation check------------------------------
        double pt_temp = 0.;
        if (isoCharged)
          for (int j = 5; j < event.size(); j++)
          {
            // only charged considered for iso cut
            if (event[j].isFinal() && event[j].isVisible() && event[j].isCharged() && j != i)
              if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        // charged + neutral pt considered for iso cut
        if (!isoCharged)
          for (int j = 5; j < event.size(); j++)
          {
            if (event[j].isFinal() && event[j].isVisible() && j != i)
              if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
                pt_temp += event[j].pT();
          }

        h->Fill(pt_temp);
        //----------------------------------------------
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_iso_All_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                                    bool isoCharged, double iso_cone_radius, double iso_pt)
{
  for (int i = 5; i < event.size(); i++)
  {
    if (event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax)
    {
      // isolation check------------------------------
      double pt_temp = 0.;
      if (isoCharged)
        for (int j = 5; j < event.size(); j++)
        {
          // only charged considered for iso cut
          if (event[j].isFinal() && event[j].isVisible() && event[j].isCharged() && j != i)
            if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
              pt_temp += event[j].pT();
        }

      // charged + neutral pt considered for iso cut
      if (!isoCharged)
        for (int j = 5; j < event.size(); j++)
        {
          if (event[j].isFinal() && event[j].isVisible() && j != i)
            if (TMath::Sqrt((event[i].phi() - event[j].phi()) * (event[i].phi() - event[j].phi()) + (event[i].eta() - event[j].eta()) * (event[i].eta() - event[j].eta())) < iso_cone_radius)
              pt_temp += event[j].pT();
        }

      h->Fill(pt_temp);
      //----------------------------------------------
    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_invXsec_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() > 90) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}


//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_TH2_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH2 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {

      int mID = event[event[event[i].iTopCopyId()].mother1()].id();

      h->Fill( electronMotherName[0], event[i].pT(), 1. );
      
      if( event[i].id() == 11 ) 
	h->Fill( electronMotherName[1], event[i].pT(), 1. );
      if( event[i].id() == -11 )
	h->Fill( electronMotherName[2], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 1000 )
	h->Fill( electronMotherName[3], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 500 &&
	  TMath::Abs(mID) < 549 )
	h->Fill( electronMotherName[4], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 400 &&
	  TMath::Abs(mID) < 439)
	h->Fill( electronMotherName[5], event[i].pT(), 1. );
      if( mID == 15 )
	h->Fill( electronMotherName[6], event[i].pT(), 1. );
      if( mID == -15 )
	h->Fill( electronMotherName[7], event[i].pT(), 1. );
      if( mID == -24 )
	h->Fill( electronMotherName[8], event[i].pT(), 1. );
      if( mID == 24 )
	h->Fill( electronMotherName[9], event[i].pT(), 1. );
      if( mID == 23 )
	h->Fill( electronMotherName[10], event[i].pT(), 1. );
      if( mID == 22 )
	h->Fill( electronMotherName[11], event[i].pT(), 1. );
      if( mID == 111 )
	h->Fill( electronMotherName[12], event[i].pT(), 1. );
      if( mID == 221 )
	h->Fill( electronMotherName[13], event[i].pT(), 1. );
      if( mID == 223 )
	h->Fill( electronMotherName[14], event[i].pT(), 1. );
      if( mID == 310 )
	h->Fill( electronMotherName[15], event[i].pT(), 1. );
      if( mID == 333 ||
	  mID == 113 ||
	  mID == 443)
	h->Fill( electronMotherName[16], event[i].pT(), 1. );

    }
  }
  return;
}
//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[event[event[i].iTopCopy()].mother1()].id() );
    }
  }
  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Fill_Electron_Pt_ByTopMotherID(Pythia8::Event &event, float etaMax, TH1 *h, std::vector <int> vec_id){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax) {
      int x = event[event[event[i].iTopCopy()].mother1()].id();
      if(std::find(vec_id.begin(), vec_id.end(), x) != vec_id.end()) {
	h->Fill( event[i].pT() );
      }
    }
  }
  return;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// void PythiaAnalysisHelper::Add_Histos_Scale_Write2File( std::vector <TH1D*>& vec, TH1* final_histo, TFile &file, double etaRange, bool useRap, bool isInvariantXsec){

//   file.cd();

//   for(unsigned int i = 0; i < vec.size(); i++){
//     final_histo->Add(vec.at(i));
//     vec.at(i)->Scale(1./etaRange, "width");
//     vec.at(i)->SetXTitle("p_{T} (GeV/#it{c})");
//     if(isInvariantXsec){
//       if(useRap) vec.at(i)->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
//       else vec.at(i)->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
//     }else{
//       if(useRap) vec.at(i)->SetYTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
//       else vec.at(i)->SetYTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
//     }
//     vec.at(i)->Write();
//   }

//   final_histo->Scale(1./etaRange, "width");
//   final_histo->SetXTitle("p_{T} (GeV/#it{c})");
//   if(isInvariantXsec){
//     if(useRap) final_histo->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
//     else final_histo->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
//   }else{
//     if(useRap) final_histo->SetYTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
//     else final_histo->SetYTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
//   }
//   final_histo->Write();


//   return;
// }

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Add_Histos_Scale_Write2File( std::vector <TH1D*>& vec, TH1* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec){

  file.cd();
  dir->cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width");
    vec.at(i)->SetXTitle("p_{T} (GeV/#it{c})");
    if(isInvariantXsec){
      if(useRap) vec.at(i)->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
      else vec.at(i)->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
    }else{
      if(useRap) vec.at(i)->SetYTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
      else vec.at(i)->SetYTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
    }

    // special treatment for gamma jet correlations histos
    TString axTitleTemp = vec.at(i)->GetName();
    if( axTitleTemp.Contains("xJetGamma") ){
      vec.at(i)->SetXTitle("x_{J#gamma}");
      vec.at(i)->SetYTitle("d#sigma_{J#gamma}/dx_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("dPhiJetGamma") ){
      vec.at(i)->SetXTitle("#Delta #phi_{J#gamma}");
      vec.at(i)->SetYTitle("d#sigma_{J#gamma}/d#Delta #phi_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("chJetTrackMult") ){
      vec.at(i)->SetXTitle("number of charged tracks per jet N_{chTr,jet}");
      vec.at(i)->SetYTitle("#sigma_{chTr,jet}/dN_{chTr,jet} (pb)");
    }
    if( axTitleTemp.Contains("xObs") ){
      vec.at(i)->SetXTitle("x^{obs}_{Pb}");
      vec.at(i)->SetYTitle("d#sigma_{J#gamma}/dx^{obs}_{Pb} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dPhi") ){
      vec.at(i)->SetXTitle("#Delta #phi_{#gamma--track}");
      vec.at(i)->SetYTitle("d#sigma_{tracks}/d#Delta #phi_{#gamma--track} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dEta") ){
      vec.at(i)->SetXTitle("#Delta #eta_{#gamma--track}");
      vec.at(i)->SetYTitle("d#sigma_{tracks}/d#Delta #eta_{#gamma--track} (pb)");
    }
    if( axTitleTemp.Contains("Bjorken") ){
      vec.at(i)->SetXTitle("real Bjorken x for given gamma-jet pair");
      vec.at(i)->SetYTitle("d#sigma_{J#gamma}/dx (pb)");
    }
    if( axTitleTemp.Contains("xSecTriggerGamma") ){
      vec.at(i)->SetXTitle("");
      vec.at(i)->SetYTitle("#sigma_{J#gamma} (pb)");
    }

    vec.at(i)->Write();

  }

  final_histo->Scale(1./etaRange, "width");
  final_histo->SetXTitle("p_{T} (GeV/#it{c})");
  if(isInvariantXsec){
    if(useRap) final_histo->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
    else final_histo->SetYTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
  }else{
    if(useRap) final_histo->SetYTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
    else final_histo->SetYTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
  }

    // special treatment for gamma jet correlations histos
    TString axTitleTemp = final_histo->GetName();
    if( axTitleTemp.Contains("xJetGamma") ){
      final_histo->SetXTitle("x_{J#gamma}");
      final_histo->SetYTitle("d#sigma_{J#gamma}/dx_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("dPhiJetGamma") ){
      final_histo->SetXTitle("#Delta #phi_{J#gamma}");
      final_histo->SetYTitle("d#sigma_{J#gamma}/d#Delta #phi_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("chJetTrackMult") ){
      final_histo->SetXTitle("number of charged tracks per jet N_{chTr,jet}");
      final_histo->SetYTitle("#sigma_{chTr,jet}/dN_{chTr,jet} (pb)");
    }
    if( axTitleTemp.Contains("xObs") ){
      final_histo->SetXTitle("x^{obs}_{Pb}");
      final_histo->SetYTitle("d#sigma_{J#gamma}/dx^{obs}_{Pb} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dPhi") ){
      final_histo->SetXTitle("#Delta #phi_{#gamma--track}");
      final_histo->SetYTitle("d#sigma_{tracks}/d#Delta #phi_{#gamma-track} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dEta") ){
      final_histo->SetXTitle("#Delta #eta_{#gamma--track}");
      final_histo->SetYTitle("d#sigma_{tracks}/d#Delta #eta_{#gamma-track} (pb)");
    }
    if( axTitleTemp.Contains("Bjorken") ){
      final_histo->SetXTitle("real Bjorken x for given gamma-jet pair");
      final_histo->SetYTitle("d#sigma_{J#gamma}/dx (pb)");
    }
    if( axTitleTemp.Contains("xSecTriggerGamma") ){
      final_histo->SetXTitle("");
      final_histo->SetYTitle("#sigma_{J#gamma} (pb)");
    }

  
  final_histo->Write();

  gROOT->cd();

  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Add_Histos_Scale_Write2File_Powheg( std::vector <TH1D>& vec, TFile &file, double invScaleFac){

  file.cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    vec.at(i).Scale(1./invScaleFac, "width");
    vec.at(i).SetXTitle("p_{T} (GeV/#it{c})");
    vec.at(i).SetYTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
    if(invScaleFac == 1.)     vec.at(i).SetYTitle("#frac{d#sigma}{dp_{T}} (pb)");

    // special treatment for gamma jet correlations histos
    TString axTitleTemp = vec.at(i).GetName();
    if( axTitleTemp.Contains("xJetGamma") ){
      vec.at(i).SetXTitle("x_{J#gamma}");
      vec.at(i).SetYTitle("d#sigma_{J#gamma}/dx_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("dPhiJetGamma") ){
      vec.at(i).SetXTitle("#Delta #phi_{J#gamma}");
      vec.at(i).SetYTitle("d#sigma_{J#gamma}/d#Delta #phi_{J#gamma} (pb)");
    }
    if( axTitleTemp.Contains("chJetTrackMult") ){
      vec.at(i).SetXTitle("number of charged tracks per jet N_{chTr,jet}");
      vec.at(i).SetYTitle("#sigma_{chTr,jet}/dN_{chTr,jet} (pb)");
    }
    if( axTitleTemp.Contains("xObs") ){
      vec.at(i).SetXTitle("x^{obs}_{Pb}");
      vec.at(i).SetYTitle("d#sigma_{J#gamma}/dx^{obs}_{Pb} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dPhi") ){
      vec.at(i).SetXTitle("#Delta #phi_{#gamma--track}");
      vec.at(i).SetYTitle("d#sigma_{tracks}/d#Delta #phi_{#gamma--track} (pb)");
    }
    if( axTitleTemp.Contains("isoCone_track_dEta") ){
      vec.at(i).SetXTitle("#Delta #eta_{#gamma--track}");
      vec.at(i).SetYTitle("d#sigma_{tracks}/d#Delta #eta_{#gamma--track} (pb)");
    }
    if( axTitleTemp.Contains("Bjorken") ){
      vec.at(i).SetXTitle("real Bjorken x for given gamma-jet pair");
      vec.at(i).SetYTitle("d#sigma_{J#gamma}/dx (pb)");
    }
    if( axTitleTemp.Contains("xSecTriggerGamma") ){
      vec.at(i).SetXTitle("");
      vec.at(i).SetYTitle("#sigma_{J#gamma} (pb)");
    }

    vec.at(i).Write();

  }

  gROOT->cd();

  return;
}

//----------------------------------------------------------------------
void PythiaAnalysisHelper::Add_Histos_Scale_Write2File( std::vector <TH2D*>& vec, TH2* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec){

  file.cd();
  dir->cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width"); // NB: width is dividing by bin area in TH2 -> here it means dividing by pt bin width, because other bin boundary is of length 1
    vec.at(i)->SetXTitle("electron Mother");
    vec.at(i)->SetYTitle("p_{T} (GeV/#it{c})");
    if(isInvariantXsec){
      if(useRap) vec.at(i)->SetZTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
      else vec.at(i)->SetZTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
    }else{
      if(useRap) vec.at(i)->SetZTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
      else vec.at(i)->SetZTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
    }
    vec.at(i)->Write();
  }

  final_histo->Scale(1./etaRange, "width");
  final_histo->SetXTitle("electron Mother");
  final_histo->SetYTitle("p_{T} (GeV/#it{c})");
  if(isInvariantXsec){
    if(useRap) final_histo->SetZTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} (pb)");
    else final_histo->SetZTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
  }else{
    if(useRap) final_histo->SetZTitle("#frac{d^{2}#sigma}{dp_{T}dy} (pb)");
    else final_histo->SetZTitle("#frac{d^{2}#sigma}{dp_{T}d#eta} (pb)");
  }
  final_histo->Write();

  gROOT->cd();

  return;
}


//----------------------------------------------------------------------
// used for powheg events
void PythiaAnalysisHelper::FillForEachWeight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights){

  if(vec_h.size() < 1){
    printf("FillForEachWeight: no histograms in vector. Aborting...\n");
    return;
  }
  for(unsigned long int i = 0; i < vec_weights.size(); i++)
    vec_h.at(i).Fill(val,vec_weights.at(i));
  
  return;
}
