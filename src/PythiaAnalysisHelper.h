#ifndef _PYTHIAANALYSISHELPER_h_included_
#define _PYTHIAANALYSISHELPER_h_included_

#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <vector>

class PythiaAnalysisHelper{

 public:

  PythiaAnalysisHelper(){};
  //  ~PythiaAnalysisHelper(){} // destructor not needed, if no member variables

  void Set_Pythia_Randomseed(Pythia8::Pythia &p); // set seed with ROOT's TRandom3
  void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv); // set stuff with one line; check .cxx for enlightenment
  void ProcessSwitch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p); // set pTHat bin specific stuff

  // fill "normal" spectra
  void Fill_Pi0_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // fill pt of all pi0 in the event within eta range
  void Fill_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // fill pt of all primary pi0 (i.e. secondary corrected)
  void Fill_Eta_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // fill pt of all eta
  void Fill_EtaPrime_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // fill pt of all eta prime
  void Fill_Omega_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // fill pt of all omega
  void Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all direct photons
  void Fill_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all shower photons
  void Fill_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all LO photons
  void Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // same with decay photons
  void Fill_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                    bool isoCharged, double iso_cone_radius, double iso_pt); // fill pt of isolated direct photons

  void Fill_Electron_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all electrons
  void Fill_ElectronNeg_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all negative electrons
  void Fill_ElectronPos_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all positive electrons = positrons

  // fill invXsec spectra
  void Fill_invXsec_Pi0_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Pi0Primary_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Eta_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_EtaPrime_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Omega_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Shower_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_222_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // same with decay photons
  void Fill_invXsec_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                             bool isoCharged, double iso_cone_radius, double iso_pt); //  ...if you want to apply 1/(pt*2pi) as weight


  // mother stuff
  void Fill_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH1 *h); // fill the top mothers id code into the histo
  void Fill_TH2_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH2 *h); // fill the top mothers id code into the histo
  void Fill_Electron_Pt_ByTopMotherID(Pythia8::Event &event, float etaMax, TH1 *h, std::vector <int> id); // fill electron by topmother id


  // post-processing
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, double etaRange, bool useRap, bool isInvariantXsec = false); //
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec = false); // adds temp. histos from pTHat bins to final; normalize for pt bin width and eta range (= no 2*Pi applied here!)
  void Add_Histos_Scale_Write2File( std::vector <TH2D*> &vec_temp_histo, TH2* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec = false); // adds temp. histos from pTHat bins to final; normalize for pt bin width and eta range (= no 2*Pi applied here!)

 private:

  const char *electronMotherName[17] = {"all",
					"neg","pos",
					"Baryons",
					"B-Mesons"    ,"D-Mesons",
					"#tau^{-}"    ,"#tau^{+}",
					"W^{-}"       ,"W^{+}",
					"Z^{0}"       ,"#gamma^{*}",
					"#pi^{0}"     ,"#eta",
					"#omega"      ,"K^{0}_{s}",
					"#Phi,#rho,J/#Psi"};

};

#endif
