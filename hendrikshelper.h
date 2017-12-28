//----------------------------------------------------------------------
// hendrikshelper.h : helpful functions for standard pythia stuff
//----------------------------------------------------------------------

#ifndef _HENDRIKSHELPER_h_included_
#define _HENDRIKSHELPER_h_included_

#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include <vector>

class HendriksHelper{

 public:

  HendriksHelper(){};
  // ~HendriksHelper(); // destructor not needed, because no member variables

  void Set_Pythia_Randomseed(Pythia8::Pythia &p); // set seed with ROOT's TRandom3
  void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv); // set stuff with one line; check .cxx for enlightenment
  void ProcessSwitch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p); // set pTHat bin specific stuff


  // fill "normal" spectra
  void Fill_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all non-decay (~direct) photons in the event within y range
  void Fill_Electron_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all electrons
  void Fill_Pi0_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all pi0
  void Fill_Eta_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all eta
  void Fill_ElectronNeg_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all negative electrons
  void Fill_ElectronPos_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // fill pt of all positive electrons = positrons
  void Fill_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                    bool isoCharged, double iso_cone_radius, double iso_pt); // fill isolated non-decay photons
  void Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // same with decay photons

  // fill invXsec spectra
  void Fill_invXsec_Direct_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Pi0_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Eta_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Direct_Iso_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h,
                                             bool isoCharged, double iso_cone_radius, double iso_pt); //  ...if you want to apply 1/(pt*2pi) as weight
  void Fill_invXsec_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // same with decay photons


  // mother stuff
  void Fill_Electron_TopMotherID(Pythia8::Event &event, float etaMax, TH1 *h); // fill the top mothers id code into the histo
  void Fill_Electron_ByTopMotherID(Pythia8::Event &event, float etaMax, TH1 *h, std::vector <int> id); // fill electron by topmother id


  // post-processing
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, double etaRange); //
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, TDirectory *dir, double etaRange); // adds temp. histos from pTHat bins to final; normalize for pt bin width and y range (= no 2*Pi applied here!)

 private:

};

#endif
