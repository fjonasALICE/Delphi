#ifndef _PYTHIAANALYSISHELPER_h_included_
#define _PYTHIAANALYSISHELPER_h_included_

#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <vector>
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "fastjet/ClusterSequence.hh"

using std::cout;
using namespace Pythia8;

using fastjet::PseudoJet;
using fastjet::JetDefinition;
using fastjet::ClusterSequence;
using fastjet::antikt_algorithm;

class PythiaAnalysisHelper{

 public:

  PythiaAnalysisHelper(){};
  //  ~PythiaAnalysisHelper(){} // destructor not needed, if no member variables

  void Set_Pythia_Randomseed(Pythia8::Pythia &p); // set seed with ROOT's TRandom3
  void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv); // set stuff with one line; check .cxx for enlightenment
  void Write_README(Pythia8::Pythia &p, TFile &file, int argc, char **argv, string pdfA, string pdfB, double *pTHatBin);
  void ProcessSwitch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p); // set pTHat bin specific stuff

  double CorrectPhiDelta(double a, double b);
  double XObs_pGoing(PseudoJet &hadronjet, PseudoJet &photonjet, double beamEnergy);
  double XObs_PbGoing(PseudoJet &hadronjet, PseudoJet &photonjet, double beamEnergy);

  bool IsPhotonIsolated(Event &event, int iPhoton, const double &etaAbsMaxPhoton, const double &isoConeRadius, const double &isoPtMax, double UEPtDensity, TH1D *h_phi, TH1D *h_eta, TH1D *h_isoPt, TH1D *h_isoPt_corrected);
  bool IsPhotonIsolatedPowheg(Event &event, int iPhoton, const double &etaAbsMaxPhoton, const double &isoConeRadius, const double &isoPtMax, double UEPtDensity, vector<TH1D> &vec_phi, vector<TH1D> &vec_eta, vector<TH1D> &vec_isoPt, vector<TH1D> &vec_isoPt_corrected, vector<double> vec_weights);

  // fill "normal" spectra
  void Fill_Electron_Pt(Pythia8::Event &event, float etaMax, bool useRap, TH1 *h);
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
  //  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, double etaRange, bool useRap, bool isInvariantXsec = false); // not used anymore, marked for deletion
  void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec = false); // adds temp. histos from pTHat bins to final; normalize for pt bin width and eta range (= no 2*Pi applied here!)
  void Add_Histos_Scale_Write2File_Powheg( std::vector <TH1D>& vec, TFile &file, double invScaleFac);
  void Add_Histos_Scale_Write2File( std::vector <TH2D*> &vec_temp_histo, TH2* final_histo, TFile &file, TDirectory *dir, double etaRange, bool useRap, bool isInvariantXsec = false); // adds temp. histos from pTHat bins to final; normalize for pt bin width and eta range (= no 2*Pi applied here!)
  void FillForEachWeight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights);


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

  static const int ptBins = 97;

  double ptBinArray[ptBins+1] = {0.0,0.2,0.4,0.6,0.8,
					1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,
					10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,
					30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,
					300};
  //  double ptBinArray[ptBins+1];
  /* for(int i=0; i < ptBins+1; i++){ */
  /*   if(i <= 4)        ptBinArray[i] = 0.2*i; // 0.0 - 0.8 */
  /*   if(i >= 5 && i <= 22)  ptBinArray[i] = 0.5*i - 1.5; // 1.0 - 10.0 */
  /*   if(i >= 23 && i <= 42) ptBinArray[i] = 1.0*i - 13.;; // 10 - 30 */
  /*   if(i >= 43)      ptBinArray[i] = 5.0*i - 185.; // 30 - 300 */
  /*   //cout << "ptBinArray[" << i << "] = " << ptBinArray[i] << endl; */
  /* } */


  const int dPhiJetGamma_nBins = 68;
  const double dPhiJetGamma_min = 0.1;
  const double dPhiJetGamma_max = 3.3;

  const int dxJetGamma_nBins = 100;
  const double dxJetGamma_min = 0.0;
  const double dxJetGamma_max = 2.5;

  const int chJetTrackMult_nBins = 50;
  const double chJetTrackMult_min = 0.5;
  const double chJetTrackMult_max = 50.5;

  const int xObs_nBins = 100;
  const double xObs_min = 0.0;
  const double xObs_max = 0.1;

  const int isoCone_track_nBins = 40;
  const double isoCone_track_min = -1;
  const double isoCone_track_max = 1.;

  
 private:

};

#endif
