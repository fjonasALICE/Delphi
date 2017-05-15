#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"

using std::cout;
using namespace Pythia8;

void Set_Pythia_Randomseed(Pythia8::Pythia &p);
void Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h);

int main(int, char **);
int main(int argc, char **argv) {

  string fileName;
  Pythia p;
  Settings& settings = p.settings;
  Pythia8::Info& info = p.info; // explicit class information to avoid confusion with TMath.....Info()

  //---read commandline args----------------------------------------
  if (argc != 6) {
    printf("Invalid number of arguments\nUsage: %s [output root file] [number of events] [cm energy in GeV] [renormScaleFac] [factorMultFac]", argv[0]);
    exit(EXIT_FAILURE);
  }

  char rootFileName[1024];
  snprintf( rootFileName, sizeof(rootFileName), "pythia_%ld_RS%.2f_FS%.2f.root", 
	    strtol(argv[3], NULL, 10), strtof(argv[4], NULL), strtof(argv[5], NULL) );
  printf("The result will be written into %s", rootFileName);
  //const char *rootFileName = argv[1]; // output file
  int nEvent = strtol(argv[2], NULL, 10); // number of events
  p.readString(Form("Beams:eCM = %f", strtof(argv[3], NULL) ));
  p.readString(Form("SigmaProcess:renormMultFac = %f", strtof(argv[4],NULL) ));
  p.readString(Form("SigmaProcess:factorMultFac = %f", strtof(argv[5],NULL) ));

  Set_Pythia_Randomseed(p);
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  
  const double ptMin = 0., ptMax = 40.;
  const int ptBins = 40;
  const double etaMax = 0.9;

  TH1D *h_shower_photons = new TH1D("h_shower_photons","shower photons", ptBins, ptMin, ptMax);
  TH1D *h_shower_photons_temp = (TH1D*)h_shower_photons->Clone("h_shower_photons_temp");

  p.readString("Next:NumberCount = 50000");

  p.readString("PartonLevel:MPI = off");
  p.readString("HadronLevel:all = off");

//TimeShower:pTminChgQ = 2.
//Check:event = off
//PartonLevel:MPI = off
//PartonLevel:ISR = off
//PartonLevel:FSR = off
//PartonLevel:Remnants = off
// PartonLevel:all = off

  // Begin event loop.
  const int pTHatBins = 10;
  double pTHatBin[pTHatBins] = {2. ,13.,16.,20.,25.,
				31.,37.,45.,53.,1000.};
 
  for (int iBin = 1; iBin < (pTHatBins-1); ++iBin) {
    settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
    settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);

    if (iBin == 0) {
      p.readString("HardQCD:all = off");
      p.readString("SoftQCD:inelastic = on");
    } else {
      p.readString("HardQCD:all = on");
      p.readString("SoftQCD:inelastic = off");
    }

    p.init();

    h_shower_photons_temp->Reset();

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate event.
      if (!p.next()) continue;
      // reject softQCD events in the hardQCD regime
      if (iBin == 0 && info.isNonDiffractive() && p.info.pTHat() > pTHatBin[1]) continue; 

      Fill_Direct_Photon_Pt(p.event, etaMax, h_shower_photons_temp);
      // End of event loop. Statistics.
    }
    p.stat();

    double sigma = p.info.sigmaGen()*1e9; // cross section (picobarn)
    double sigma_per_event = sigma/p.info.weightSum();

    h_shower_photons->Add(h_shower_photons_temp, sigma_per_event);
  }

  h_shower_photons->Scale(1.,"width");

  // write to root file
  TFile file(rootFileName, "RECREATE");
  h_shower_photons->Write();
  file.Close();

  return 0;
}

void Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill(event[i].pT());
    }
  }
}

void Set_Pythia_Randomseed(Pythia8::Pythia &p){

  TRandom3 *rand = new TRandom3(0);
  p.readString("Random:setSeed = on");
  p.readString(Form("Random:seed = %ld", (long)(rand->Rndm()*9e8) ));
  delete rand;
  return;
}
