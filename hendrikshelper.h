//----------------------------------------------------------------------
// hendriks_helper.h : helpful functions for standard pythia stuff
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
  void Fill_Non_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // returns pt of all non-decay (~direct) photons in the event within eta range
  void Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h); // same with decay photons
  void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv); // set stuff with one line; check .cxx for enlightenment
  void SoftQCD_HardQCD_Switch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p, int &nEvent); // set pTHat bin specific stuff
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, double etaRange); // adds temp. histos from pTHat bins to final; normalize

 private:
  
};

#endif
