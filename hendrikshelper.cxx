//----------------------------------------------------------------------
// hendriks_helper.cxx
//----------------------------------------------------------------------
#include "hendrikshelper.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void HendriksHelper::Set_Pythia_Randomseed(Pythia8::Pythia &p){

  TRandom3 *rand = new TRandom3(0);
  p.readString("Random:setSeed = on");
  p.readString(Form("Random:seed = %ld", (long)(rand->Rndm()*9e8) ));
  delete rand;
  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv){

  p.readString(Form("Beams:eCM = %f", strtof(argv[4], NULL) ));
  if(argc > 6)
    p.readString(Form("SigmaProcess:renormMultFac = %f", strtof(argv[5],NULL) ));
  if(argc > 7)
    p.readString(Form("SigmaProcess:factorMultFac = %f", strtof(argv[6],NULL) ));

  if(argc > 5){
    if( !strcmp(argv[5],"noHadro") ){
      p.readString("HadronLevel:all = off");
      printf("---pythia events are generated without Hadronization---\n");
    }
    else if( !strcmp(argv[5],"noMPI") ){
      p.readString("PartonLevel:MPI = off");
      printf("---pythia events are generated without Multiparton Interaction---\n");
    }
    else if( !strcmp(argv[5],"noMPInoHadro") ){
      p.readString("PartonLevel:MPI = off");
      p.readString("HadronLevel:all = off");
      printf("---pythia events are generated without Multiparton Interaction and Hadronization---\n");
    }
    else if( !strcmp(argv[5],"noShower") ){
      p.readString("Check:event = off");
      p.readString("PartonLevel:MPI = off");
      p.readString("PartonLevel:FSR = off");
      p.readString("PartonLevel:ISR = off");
      p.readString("PartonLevel:Remnants = off");
      p.readString("HadronLevel:all = off");
      printf("---pythia events are generated without parton shower (only LO scattering)---\n");
    }
    else if( !strcmp(argv[5],"fullEvents") ){
      printf("---full pythia events will be generated (default)---\n");
    }
    else
      printf("---no sensible argument argv[5] is given -> full pythia events will be generated (default)---\n");
  }
  
  return;
}

//----------------------------------------------------------------------
void HendriksHelper::ProcessSwitch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p){

  if( !strcmp(argv[2],"MB") || !strcmp(argv[2],"MBVeto")){
    p.readString("SoftQCD:inelastic= on");
    return;
  }
  else if( !strcmp(argv[2],"JJ") ){
    p.readString("HardQCD:all = on");
    p.readString("HardQCD:hardccbar = on");
    p.readString("HardQCD:hardbbbar = on");
  }
  else if( !strcmp(argv[2],"PromptPhoton") ){
    p.readString("PromptPhoton:all = on");
  }
  else if( !strcmp(argv[2],"WeakBoson") ){
    p.readString("WeakBosonExchange:all = on");
    p.readString("WeakSingleBoson:all = on");
    p.readString("WeakDoubleBoson:all = on");
  }
  else
    printf("No process switched on. Provide \"MB\" or \"MBVeto\" or \"JJ\" or \"PromptPhoton\" or \"WeakBoson\" as second argument");
  
  p.settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
  p.settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);

  return;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void HendriksHelper::Fill_Pi0_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
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
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Eta_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 221 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Omega_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 223 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90)
	if( TMath::Abs(event[event[i].iTopCopy()].status() ) > 40 )
	  h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
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
void HendriksHelper::Fill_Electron_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_ElectronNeg_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_ElectronPos_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == -11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill(event[i].pT());
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
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
void HendriksHelper::Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
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
void HendriksHelper::Fill_invXsec_Pi0_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 111 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_invXsec_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
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
      h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
	}
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_invXsec_Eta_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 221 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_invXsec_Omega_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].id() == 223 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_invXsec_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_invXsec_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
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
void HendriksHelper::Fill_invXsec_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
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
void HendriksHelper::Fill_invXsec_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
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
void HendriksHelper::Fill_invXsec_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() > 90) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}


//----------------------------------------------------------------------
void HendriksHelper::Fill_TH2_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH2 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {

      int mID = event[event[event[i].iTopCopyId()].mother1()].id();
      //      std::cout << "MotherID(iTopCopyId): " << mID << std::endl;

      h->Fill( electronMotherName[0], event[i].pT(), 1. );
      
      if( event[i].id() == 11 ) 
	h->Fill( electronMotherName[1], event[i].pT(), 1. );
      if( event[i].id() == -11 )
	h->Fill( electronMotherName[2], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 1000 )
	h->Fill( electronMotherName[3], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 500 &&
	  TMath::Abs(mID) < 600 )
	h->Fill( electronMotherName[4], event[i].pT(), 1. );
      if( TMath::Abs(mID) > 400 &&
	  TMath::Abs(mID) < 500 )
	h->Fill( electronMotherName[5], event[i].pT(), 1. );
      if( mID == 11 ){
	// event.list();
	// std::cout << "electron at index " << i << std::endl;
	// std::cout << "iTopCopy: " << event[i].iTopCopy() << std::endl;
	// std::cout << "iTopCopyId: " << event[i].iTopCopyId() << std::endl;
	h->Fill( electronMotherName[6], event[i].pT(), 1. );
      }
      if( mID == -11 ){
	// event.list();
	// std::cout << "positron at index " << i << std::endl;
	// std::cout << "iTopCopy: " << event[i].iTopCopy() << std::endl;
	// std::cout << "iTopCopyId: " << event[i].iTopCopyId() << std::endl;
	h->Fill( electronMotherName[7], event[i].pT(), 1. );
      }
      if( mID == 15 )
	h->Fill( electronMotherName[8], event[i].pT(), 1. );
      if( mID == -15 )
	h->Fill( electronMotherName[9], event[i].pT(), 1. );
      if( mID == -24 )
	h->Fill( electronMotherName[10], event[i].pT(), 1. );
      if( mID == 24 )
	h->Fill( electronMotherName[11], event[i].pT(), 1. );
      if( mID == 23 )
	h->Fill( electronMotherName[12], event[i].pT(), 1. );
      if( mID == 22 )
	h->Fill( electronMotherName[13], event[i].pT(), 1. );
      if( mID == 111 )
	h->Fill( electronMotherName[14], event[i].pT(), 1. );
      if( mID == 221 )
	h->Fill( electronMotherName[15], event[i].pT(), 1. );
      if( mID == 223 )
	h->Fill( electronMotherName[16], event[i].pT(), 1. );
      if( mID == 310 )
	h->Fill( electronMotherName[17], event[i].pT(), 1. );
      if( mID == 333 ||
	  mID == 113 ||
	  mID == 443)
	h->Fill( electronMotherName[18], event[i].pT(), 1. );
      
    }
  }
  return;
}
//----------------------------------------------------------------------
void HendriksHelper::Fill_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[event[event[i].iTopCopy()].mother1()].id() );
    }
  }
  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Fill_Electron_Pt_ByTopMotherID(Pythia8::Event &event, float etaMax, TH1 *h, std::vector <int> vec_id){
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
void HendriksHelper::Add_Histos_Scale_Write2File( std::vector <TH1D*>& vec, TH1* final_histo, TFile &file, double etaRange){

  file.cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width");
    vec.at(i)->SetXTitle("p_{T} (GeV/#it{c})");
    vec.at(i)->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
    vec.at(i)->Write();
  }

  final_histo->SetXTitle("p_{T} (GeV/#it{c})");
  final_histo->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
  final_histo->Scale(1./etaRange, "width");
  final_histo->Write();


  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Add_Histos_Scale_Write2File( std::vector <TH1D*>& vec, TH1* final_histo, TFile &file, TDirectory *dir, double etaRange){

  file.cd();
  dir->cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width");
    vec.at(i)->SetXTitle("p_{T} (GeV/#it{c})");
    vec.at(i)->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
    vec.at(i)->Write();
  }

  final_histo->SetXTitle("p_{T} (GeV/#it{c})");
  final_histo->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
  final_histo->Scale(1./etaRange, "width");
  final_histo->Write();

  gROOT->cd();

  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Add_Histos_Scale_Write2File( std::vector <TH2D*>& vec, TH2* final_histo, TFile &file, TDirectory *dir, double etaRange){

  file.cd();
  dir->cd();

  for(unsigned int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width");
    vec.at(i)->SetXTitle("electron Mother");
    vec.at(i)->SetYTitle("p_{T} (GeV/#it{c})");
    vec.at(i)->SetZTitle("#frac{d#sigma}{dp_{T}d#eta}");
    vec.at(i)->Write();
  }

  final_histo->SetXTitle("electron Mother");
  final_histo->SetYTitle("p_{T} (GeV/#it{c})");
  final_histo->SetZTitle("#frac{d#sigma}{dp_{T}d#eta}");
  final_histo->Write();

  gROOT->cd();

  return;
}
