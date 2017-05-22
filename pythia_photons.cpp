//#include <iostream>
//#include <cmath>
//#include <cstring>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
//#include "TMath.h"
//#include "TH1.h"
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


  // ---for direct photons---
  // const int pTHatBins = 8;
  // double pTHatBin[pTHatBins+1] = {2., 5., 8., 12., 17.,
  //  				  23., 30., 37., 10000.};
  
  // ---for shower photons---
  const int pTHatBins = 7;
  printf("-----------------------\nusing %d pTHat bins\n-----------------------", pTHatBins);
  
  // test how much hardQCD overestimates cross section
  double pTHatBin[pTHatBins+1] = {2., 5., 8., 12., 17., 23., 30., 1000.};
  // test #2
  //double pTHatBin[pTHatBins+1] = {2., 8., 12., 17., 23., 30., 1000.};

  TH1D *h_non_decay_photons_etaTPC = new TH1D("h_non_decay_photons_etaTPC","non_decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_non_decay_photons_etaEMCal = new TH1D("h_non_decay_photons_etaEMCal","non_decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_non_decay_photons_etaPHOS = new TH1D("h_non_decay_photons_etaPHOS","non_decay photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_decay_photons_etaTPC = new TH1D("h_decay_photons_etaTPC","decay photons_etaTPC", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaEMCal = new TH1D("h_decay_photons_etaEMCal","decay photons_etaEMCal", ptBins, ptMin, ptMax);
  TH1D *h_decay_photons_etaPHOS = new TH1D("h_decay_photons_etaPHOS","decay photons_etaPHOS", ptBins, ptMin, ptMax);

  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", ptBins, ptMin, ptMax);

  // organise pTHat wise histograms in vectors
  vector <TH1D*> vec_non_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_non_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_non_decay_photons_etaPHOS_bin;

  vector <TH1D*> vec_decay_photons_etaTPC_bin;
  vector <TH1D*> vec_decay_photons_etaEMCal_bin;
  vector <TH1D*> vec_decay_photons_etaPHOS_bin;

  vector <TH1D*> vec_pTHat_bin;


  for(int i = 0; i < pTHatBins; i++){
    char buffer[64];

    int n = sprintf(buffer, "h_non_decay_photons_etaTPC_bin_%d", i) ;
    vec_non_decay_photons_etaTPC_bin.push_back( (TH1D*)h_non_decay_photons_etaTPC->Clone(buffer) );

    n = sprintf(buffer, "h_non_decay_photons_etaEMCal_bin_%d", i) ;
    vec_non_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_non_decay_photons_etaEMCal->Clone(buffer) );

    n = sprintf(buffer, "h_non_decay_photons_etaPHOS_bin_%d", i) ;
    vec_non_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_non_decay_photons_etaPHOS->Clone(buffer) );

    n = sprintf(buffer, "h_decay_photons_etaTPC_bin_%d", i) ;
    vec_decay_photons_etaTPC_bin.push_back( (TH1D*)h_decay_photons_etaTPC->Clone(buffer) );

    n = sprintf(buffer, "h_decay_photons_etaEMCal_bin_%d", i) ;
    vec_decay_photons_etaEMCal_bin.push_back( (TH1D*)h_decay_photons_etaEMCal->Clone(buffer) );

    n = sprintf(buffer, "h_decay_photons_etaPHOS_bin_%d", i) ;
    vec_decay_photons_etaPHOS_bin.push_back( (TH1D*)h_decay_photons_etaPHOS->Clone(buffer) );

    n = sprintf(buffer, "h_pTHat_bin_%d", i) ;
    vec_pTHat_bin.push_back( (TH1D*)h_pTHat->Clone(buffer) );

  }

  //--- begin pTHat bin loop ----------------------------------
  for (int iBin = 0; iBin < pTHatBins; ++iBin) {

    pyHelp.SoftQCD_HardQCD_Switch(iBin, pTHatBin, argv, p, nEvent);
    p.init();

    //--- begin event loop ----------------------------------------------
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      // Generate event.
      if (!p.next()) continue;
      // reject non-diffractive softQCD events in the hardQCD regime
      if (iBin == 0 && p.info.isNonDiffractive() && p.info.pTHat() > pTHatBin[1]) continue;
      
      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaTPC, vec_non_decay_photons_etaTPC_bin.at(iBin));
      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaEMCal, vec_non_decay_photons_etaEMCal_bin.at(iBin));
      pyHelp.Fill_Non_Decay_Photon_Pt(p.event, etaPHOS, vec_non_decay_photons_etaPHOS_bin.at(iBin));

      if( !strcmp(argv[2],"decay") ){
	pyHelp.Fill_Decay_Photon_Pt(p.event, etaTPC, vec_decay_photons_etaTPC_bin.at(iBin));
	pyHelp.Fill_Decay_Photon_Pt(p.event, etaEMCal, vec_decay_photons_etaEMCal_bin.at(iBin));
	pyHelp.Fill_Decay_Photon_Pt(p.event, etaPHOS, vec_decay_photons_etaPHOS_bin.at(iBin));
      }
      
      vec_pTHat_bin.at(iBin)->Fill(p.info.pTHat());

    }// end of event loop

    p.stat();

    double sigma = p.info.sigmaGen()*1e9; // cross section in picobarn
    double sigma_per_event = sigma/p.info.weightSum(); // weightSum = number of events in standard Pythia8
    vec_non_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma_per_event);
    vec_non_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma_per_event);
    vec_non_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma_per_event);
    vec_decay_photons_etaTPC_bin.at(iBin)->Scale(sigma_per_event);
    vec_decay_photons_etaEMCal_bin.at(iBin)->Scale(sigma_per_event);
    vec_decay_photons_etaPHOS_bin.at(iBin)->Scale(sigma_per_event);
    vec_pTHat_bin.at(iBin)->Scale(sigma_per_event);

  }// end of pTHat bin loop

  //--- write to root file ---------------------------------------
  TFile file(rootFileName, "RECREATE");

  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaTPC_bin, h_non_decay_photons_etaTPC, file, 2*etaTPC);
  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaEMCal_bin, h_non_decay_photons_etaEMCal, file, 2*etaEMCal);
  pyHelp.Add_Histos_Scale_Write2File( vec_non_decay_photons_etaPHOS_bin, h_non_decay_photons_etaPHOS, file, 2*etaPHOS);

  if( !strcmp(argv[2],"decay") ){
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaTPC_bin, h_decay_photons_etaTPC, file, 2*etaTPC);
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaEMCal_bin, h_decay_photons_etaEMCal, file, 2*etaEMCal);
    pyHelp.Add_Histos_Scale_Write2File( vec_decay_photons_etaPHOS_bin, h_decay_photons_etaPHOS, file, 2*etaPHOS);
  }
  
  pyHelp.Add_Histos_Scale_Write2File( vec_pTHat_bin, h_pTHat, file, 1.);

  file.Close();

  return 0;
}
