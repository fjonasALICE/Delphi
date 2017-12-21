#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "hendrikshelper.h"
#include "TTree.h"

using std::cout;
using namespace Pythia8;

int main(int, char **);
int main(int argc, char **argv) {

  //--- read commandline args ----------------------------------------
  if (argc < 2) {
    printf("Invalid number of arguments, needs two params: string and number of events");
    exit(EXIT_FAILURE);
  }

  Pythia p;
  HendriksHelper pyHelp;

  const int startBin = 0;

  //--- define output root file --------------------------------------
  char rootFileName[1024];
  snprintf( rootFileName, sizeof(rootFileName), "pythia_electrons_%s_%s.root", argv[1], argv[2]);
  printf("------------------------------------------------------------------------------------------\n");
  printf("The result will be written into %s\n", rootFileName);
  printf("------------------------------------------------------------------------------------------\n");

  int nEvent = strtol(argv[2], NULL, 10); // number of events

  //  pyHelp.Pass_Parameters_To_Pythia(p, argc, argv); // which energy, scales, optional master switches

  p.readString("Next:NumberCount = 100000");
  pyHelp.Set_Pythia_Randomseed(p);

  p.readString("Beams:eCM = 8000");

  p.readString("PartonLevel:MPI = off");

  p.readString("HardQCD:all = on");

  p.readString("HardQCD:hardccbar = on");
  p.readString("HardQCD:hardbbbar = on");

  p.readString("PromptPhoton:all = on");

  p.readString("WeakBosonExchange:all = on");

  p.readString("WeakSingleBoson:all = on");
  p.readString("WeakDoubleBoson:all = on");



  //--- Histograms ---------------------------------------------------
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const double ptMin = 0., ptMax = 300.;
  const int ptBins = 300;
  const double etaMax = 0.9;

  const int pTHatBins = 19;
  double pTHatBin[pTHatBins+1] = { 0.,
  				   9.  , 12. , 16. , 21. , 28.,
     				   36. , 45. , 57. , 70. , 85.,
  				   99. , 115., 132., 150., 169.,
  				   190., 212., 235 ,
				   10000. };
  // const int pTHatBins = 5;
  // double pTHatBin[pTHatBins+1] = { 20., 50., 100., 200., 400., 1000. }; 

  printf("-----------------------\nusing %d pTHat bins\n-----------------------\n", pTHatBins);

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------
  //

  TH1I *h_processCodes = new TH1I("h_processCodes", "process codes", 160, 99.5, 259.5);
  TH1I *h_electron_TopMotherID = new TH1I("h_electron_TopMotherID", "top mother of the electron", 12001, -6000.5, 6000.5);

  TH1D *h_electron       = new TH1D("h_electron","inclusive electrons", ptBins, ptMin, ptMax);
  TH1D *h_electron_QCD   = new TH1D("h_electron_QCD","electrons from HardQCD", ptBins, ptMin, ptMax);
  TH1D *h_electron_Gamma = new TH1D("h_electron_Gamma","electrons from PromptPhoton processes", ptBins, ptMin, ptMax);
  TH1D *h_electron_HF    = new TH1D("h_electron_HF","electrons from prompt HF processes", ptBins, ptMin, ptMax);
  TH1D *h_electron_Weak  = new TH1D("h_electron_Weak","electrons from prompt Weak bosons", ptBins, ptMin, ptMax);

  TH1D *h_electron_Neg   = new TH1D("h_electron_Neg","inclusive electrons (pos)", ptBins, ptMin, ptMax);
  TH1D *h_electron_Pos   = new TH1D("h_electron_Pos","inclusive electrons (neg)", ptBins, ptMin, ptMax);

  TH1D *h_electron_Z     = new TH1D("h_electron_Z","electrons from Z^{0}", ptBins, ptMin, ptMax);
  TH1D *h_electron_Wp    = new TH1D("h_electron_Wp","inclusive electrons from W^{+}", ptBins, ptMin, ptMax);
  TH1D *h_electron_Wm    = new TH1D("h_electron_Wm","inclusive electrons from W^{-}", ptBins, ptMin, ptMax);



  //----------------------------------------------------------------------------------------------------
  // check underlying born kt to see, e.g., if HardQCD cross section does not blow up
  TH1D *h_pTHat = new TH1D("h_pTHat","pTHat aka born kt", 1000, ptMin, ptMax);
  // store the weight sum for proper normalization afterwards
  TH1D *h_weightSum = new TH1D("h_weightSum","sum of weights", 1, 0., 1.);

  // organise pTHat wise histograms in vectors
  vector <TH1D*> vec_electron_bin,
    vec_electron_QCD_bin,
    vec_electron_Gamma_bin,
    vec_electron_HF_bin,
    vec_electron_Weak_bin,
    vec_electron_Neg_bin,
    vec_electron_Pos_bin,
    vec_electron_Z_bin,
    vec_electron_Wp_bin,
    vec_electron_Wm_bin;

  vector <TH1I*> vec_processCodes_bin;
  vector <TH1I*> vec_electron_TopMotherID_bin;

  //----------------------------------------------------------------------------------------------------

  vector <TH1D*> vec_pTHat_bin;
  vector <TH1D*> vec_weightSum_bin;


  for(int i = 0; i < pTHatBins; i++){

    // electron histos
    vec_electron_bin.push_back( (TH1D*)h_electron->Clone(Form("h_electron_bin_%02d",i)) );
    vec_electron_QCD_bin.push_back( (TH1D*)h_electron_QCD->Clone(Form("h_electron_QCD_bin_%02d",i)) );
    vec_electron_Gamma_bin.push_back( (TH1D*)h_electron_Gamma->Clone(Form("h_electron_Gamma_bin_%02d",i)) );
    vec_electron_HF_bin.push_back( (TH1D*)h_electron_HF->Clone(Form("h_electron_HF_bin_%02d",i)) );
    vec_electron_Weak_bin.push_back( (TH1D*)h_electron_Weak->Clone(Form("h_electron_Weak_bin_%02d",i)) );
    vec_electron_Neg_bin.push_back( (TH1D*)h_electron_Neg->Clone(Form("h_electron_Neg_bin_%02d",i)) );
    vec_electron_Pos_bin.push_back( (TH1D*)h_electron_Pos->Clone(Form("h_electron_Pos_bin_%02d",i)) );
    vec_electron_Z_bin.push_back( (TH1D*)h_electron_Z->Clone(Form("h_electron_Z_bin_%02d",i)) );
    vec_electron_Wp_bin.push_back( (TH1D*)h_electron_Wp->Clone(Form("h_electron_Wp_bin_%02d",i)) );
    vec_electron_Wm_bin.push_back( (TH1D*)h_electron_Wm->Clone(Form("h_electron_Wm_bin_%02d",i)) );

    vec_processCodes_bin.push_back( (TH1I*)h_processCodes->Clone(Form("h_processCodes_bin_%02d",i)) );
    vec_electron_TopMotherID_bin.push_back( (TH1I*)h_electron_TopMotherID->Clone(Form("h_electron_TopMotherID_bin_%02d",i)) );

    //----------------------------------------------------------------------------------------------------
    vec_weightSum_bin.push_back( (TH1D*)h_weightSum->Clone(Form( "h_weightSum_bin_%02d", i )) );
    vec_pTHat_bin.push_back( (TH1D*)h_pTHat->Clone(Form("h_pTHat_bin_%02d",i)) );

  }

  //--- begin pTHat bin loop ----------------------------------
  for (int iBin = startBin; iBin < pTHatBins; ++iBin) {

    p.settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
    p.settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);
    //    pyHelp.SoftQCD_HardQCD_Switch(iBin, pTHatBin, argv, p, nEvent);
    p.init();

    //--- begin event loop ----------------------------------------------
    for (int iEvent = 1; iEvent <= nEvent; ++iEvent) {
      // Generate event.
      if (!p.next()) continue;

      if(iEvent == 1)
        cout << "energy of beam a = " << p.event[1].e() << endl
             << "energy of beam b = " << p.event[2].e() << endl;


      pyHelp.Fill_Electron_Pt(p.event, etaMax, vec_electron_bin.at(iBin));
      pyHelp.Fill_ElectronNeg_Pt(p.event, etaMax, vec_electron_Neg_bin.at(iBin));
      pyHelp.Fill_ElectronPos_Pt(p.event, etaMax, vec_electron_Pos_bin.at(iBin));

      std::vector<int> id_Z{23};
      std::vector<int> id_Wp{24};
      std::vector<int> id_Wm{-24};

      pyHelp.Fill_Electron_ByTopMotherID(p.event, etaMax, vec_electron_Z_bin.at(iBin), id_Z);
      pyHelp.Fill_Electron_ByTopMotherID(p.event, etaMax, vec_electron_Wp_bin.at(iBin), id_Wp);
      pyHelp.Fill_Electron_ByTopMotherID(p.event, etaMax, vec_electron_Wm_bin.at(iBin), id_Wm);


      if( p.info.code() >= 111 && p.info.code() <= 116)
	pyHelp.Fill_Electron_Pt(p.event, etaMax, vec_electron_QCD_bin.at(iBin));

      if( p.info.code() >= 201 && p.info.code() <= 205)
	pyHelp.Fill_Electron_Pt(p.event, etaMax, vec_electron_Gamma_bin.at(iBin));

      if( p.info.code() >= 121 && p.info.code() <= 124)
	pyHelp.Fill_Electron_Pt(p.event, etaMax, vec_electron_HF_bin.at(iBin));

      if( p.info.code() >= 211 && p.info.code() <= 254)
	pyHelp.Fill_Electron_Pt(p.event, etaMax, vec_electron_Weak_bin.at(iBin));

      vec_processCodes_bin.at(iBin)->Fill(p.info.code());
      pyHelp.Fill_Electron_TopMotherID(p.event, etaMax, vec_electron_TopMotherID_bin.at(iBin));

      //----------------------------------------------------------------------------------------------------
      vec_pTHat_bin.at(iBin)->Fill(p.info.pTHat());

    }// end of event loop

    p.stat();

    double sigma = p.info.sigmaGen()*1e9; // cross section in picobarn
    //    double sigma_per_event = sigma/p.info.weightSum(); // weightSum = number of events in standard Pythia8

    vec_weightSum_bin.at(iBin)->SetBinContent(1,p.info.weightSum());
    h_weightSum->SetBinContent(1,h_weightSum->GetBinContent(1)+p.info.weightSum());
    cout << "- - - weightSum() = " << p.info.weightSum() << endl;

    vec_electron_bin.at(iBin)->Scale(sigma);
    vec_electron_QCD_bin.at(iBin)->Scale(sigma);
    vec_electron_Gamma_bin.at(iBin)->Scale(sigma);
    vec_electron_HF_bin.at(iBin)->Scale(sigma);
    vec_electron_Weak_bin.at(iBin)->Scale(sigma);

    vec_electron_Neg_bin.at(iBin)->Scale(sigma);
    vec_electron_Pos_bin.at(iBin)->Scale(sigma);
    vec_electron_Z_bin.at(iBin)->Scale(sigma);
    vec_electron_Wp_bin.at(iBin)->Scale(sigma);
    vec_electron_Wm_bin.at(iBin)->Scale(sigma);

    //----------------------------------------------------------------------------------------------------
    vec_pTHat_bin.at(iBin)->Scale(sigma);


  }// end of pTHat bin loop

  //--- write to root file ---------------------------------------
  TFile file(rootFileName, "RECREATE");


  pyHelp.Add_Histos_Scale_Write2File( vec_electron_bin,       h_electron,       file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_QCD_bin,   h_electron_QCD,   file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Gamma_bin, h_electron_Gamma, file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_HF_bin,    h_electron_HF,    file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Weak_bin,  h_electron_Weak,  file, 2*etaMax);

  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Neg_bin,   h_electron_Neg,   file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Pos_bin,   h_electron_Pos,   file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Z_bin,     h_electron_Z,     file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Wp_bin,    h_electron_Wp,    file, 2*etaMax);
  pyHelp.Add_Histos_Scale_Write2File( vec_electron_Wm_bin,    h_electron_Wm,    file, 2*etaMax);

  //----------------------------------------------------------------------------------------------------
  pyHelp.Add_Histos_Scale_Write2File( vec_pTHat_bin, h_pTHat, file, 1.);

  for(int iBin=0; iBin < pTHatBins; iBin++){
    vec_weightSum_bin.at(iBin)->Write();
    vec_processCodes_bin.at(iBin)->Write();
    vec_electron_TopMotherID_bin.at(iBin)->Write();
  }
  h_weightSum->Write();

  file.Close();

  return 0;
}
