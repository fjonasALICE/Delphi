//----------------------------------------------------------------------
// hendriks_helper.cxx
//----------------------------------------------------------------------
#include "hendrikshelper.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"


//----------------------------------------------------------------------
void HendriksHelper::Fill_Non_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill(event[i].pT());
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
void HendriksHelper::Fill_invXsec_Non_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill( event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Fill_Non_Decay_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
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
void HendriksHelper::Fill_invXsec_Non_Decay_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
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
void HendriksHelper::Fill_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax ) {
      h->Fill( event[event[event[i].iTopCopy()].mother1()].id() );
    }
  }
  return;
}

//----------------------------------------------------------------------
void HendriksHelper::Fill_Electron_ByTopMotherID(Pythia8::Event &event, float etaMax, TH1 *h, std::vector <int> vec_id){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && TMath::Abs(event[i].id()) == 11 && TMath::Abs(event[i].eta()) < etaMax) {
      int x = event[event[event[i].iTopCopy()].mother1()].id();
      if(std::find(vec_id.begin(), vec_id.end(), x) != vec_id.end()) {
	std::cout << "wuff" << std::endl;
	h->Fill( event[i].pT() );
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
void HendriksHelper::Fill_invXsec_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() > 90) h->Fill(event[i].pT(), 1./event[i].pT()/(2*TMath::Pi()) );
    }
  }
  return;
}

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
  p.readString(Form("SigmaProcess:renormMultFac = %f", strtof(argv[5],NULL) ));
  p.readString(Form("SigmaProcess:factorMultFac = %f", strtof(argv[6],NULL) ));

  if(argc > 7){
    if( !strcmp(argv[7],"noHadro") ){
      p.readString("HadronLevel:all = off");
    }
    else if( !strcmp(argv[7],"noMPI") ){
      p.readString("PartonLevel:MPI = off");
    }
    else if( !strcmp(argv[7],"noMPInoHadro") ){
      p.readString("PartonLevel:MPI = off");
      p.readString("HadronLevel:all = off");
    }
    else if( !strcmp(argv[7],"noShower") ){
      p.readString("Check:event = off");
      p.readString("PartonLevel:MPI = off");
      p.readString("PartonLevel:FSR = off");
      p.readString("PartonLevel:ISR = off");
      p.readString("PartonLevel:Remnants = off");
      p.readString("HadronLevel:all = off");
    }
  }
  else printf("---no optional argument is given -> full pythia events will be generated---\n");

  return;
}

//----------------------------------------------------------------------
void HendriksHelper::SoftQCD_HardQCD_Switch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p, int &nEvent){

  if( !strcmp(argv[2],"frag") ||
      !strcmp(argv[2],"decay") ){

    // SoftQCD only in the first pthat bin, HardQCD otherwise
    if (iBin == 0) {
      p.readString("HardQCD:all = off");
      p.readString("SoftQCD:inelastic = on");
      nEvent *= 50;
    } else {
      if(iBin == 1)
	//nEvent /= 50;
      p.readString("HardQCD:all = on");
      p.readString("SoftQCD:inelastic = off");
    }

  }
  else if( !strcmp(argv[2],"dir") ){
    p.readString("PromptPhoton:qg2qgamma = on");
    p.readString("PromptPhoton:qqbar2ggamma = on");
  }
  else
    printf("No process switched on. Provide \"dir\" or \"frag\" or \"decay\" as second argument");

  p.settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
  p.settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);


  return;
}

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
