//----------------------------------------------------------------------
// example code of pythia8 shower applied on lhef files
// for the powheg directphoton process
//----------------------------------------------------------------------
// There are analysis cuts applied following an ATLAS measurement
// of isolated photons (INSPIREHEP 1510441).
// The output should reproduce fig.2 from powheg direct photon
// simulations from INSPIREHEP 1623205.
//----------------------------------------------------------------------
// author: H. Poppenborg (hendrik.poppenborg@wwu.de)
//----------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <map>
#include <cstring>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "ShowerAnalysis.h"
#include "PythiaAnalysisHelper.h"

#include "fastjet/ClusterSequence.hh"

#include "QEDQCDPowhegHooks.h"

void ReadInWeightIDs(Pythia &p, vector<string> &vec_weightsID, bool &isSudaWeight);
double GetUEPtDensity(Event &event, int iPhoton);

int main(int, char **);
int main(int argc, char **argv) {

  string fileName; Pythia p; bool loadhooks;
  QEDQCDPowhegHooks *powhegHooks = 0; // POWHEG UserHooks
  p.readFile("shower.conf");

  //---read commandline args----------------------------------------
  if (argc < 3) {
    cout << endl << "Usage: " << argv[0]
         << "outputfile.root boostFac eventfile.lhe" << endl;
    exit(EXIT_FAILURE);
  } 
  const char *rootFileName = argv[1]; // output file

  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  vector<PseudoJet> vJets;
  ClusterSequence *cs = 0;
  
  const int nPtBins = 1000;
  const double PtBinsMin = 0., PtBinsMax = 100.;
  // vectors of histograms for different weights (e.g. for scale/pdf variation)
  vector<TH1D> vec_directphoton_pt; // = no decay photons
  vector<TH1D> vec_directphoton_pt_leading; // only hardest direct photon in event
  vector<TH1D> vec_directphoton_pt_leading_bornveto00; // born veto off
  vector<TH1D> vec_directphoton_pt_leading_bornveto20; // hard born veto
  vector<TH1D> vec_directphoton_pt_leading_bornveto40; // soft born veto
  vector<TH1D> vec_isodirectphoton_pt, vec_isodirectphoton_pt_leading;
  vector<TH1D> vec_chjet_pt, vec_chjet_pt_leading;
  vector<TH1D> vec_isoCone_track_phi, vec_isoCone_track_eta;
  vector<TH1D> vec_dPhiJetGamma, vec_dPhiJetGamma_noDeltaPhiCut;
  vector<TH1D> vec_xJetGamma;
  vector<TH1D> vec_chJetTrackMult;
  vector<TH1D> vec_xObs_pGoing, vec_xObs_PbGoing;
  vector<TH1D> vec_xBjorken_1, vec_xBjorken_2;
  vector<TH1D> vec_xBjorken_1_PDF, vec_xBjorken_2_PDF;
  vector<TH1D> vec_xSecTriggerGamma;
  vector<TH1D> vec_isoPt, vec_UEPtDensity, vec_isoPt_corrected;
  TH1D *h_nEvents = new TH1D("h_nEvents", "number of events", 4, 0.5, 4.5);
  h_nEvents->GetXaxis()->SetBinLabel(1, "bornveto2.5(std)");
  h_nEvents->GetXaxis()->SetBinLabel(2, "bornveto0.0(off)");
  h_nEvents->GetXaxis()->SetBinLabel(3, "bornveto2.0");
  h_nEvents->GetXaxis()->SetBinLabel(4, "bornveto3.0");
  
  // prepare bookkeeping of weights
  //----------------------------------------------------------------------
  bool isSudaWeight = false; // was photon radiation enhanced?
  double sudaWeight = 1.;    // reweighting factor associated with radiation enhancement
  vector<double> vec_weights;   // shall later contain: sudaWeight * primary event weight (using vector to store multiple weights, e.g for scale/pdf variation)
  vector<string> vec_weightsID;// vector storing descriptive id of weights

  // pythia settings required for usage with powheg
  //----------------------------------------------------------------------
  p.readString("Beams:frameType = 4");
  p.readString("Next:numberCount = 50000");

  double boostBetaZ = strtof(argv[2], NULL);
  if(boostBetaZ != 0.) printf("\nApplying a boost along z direction with beta_z = %f\n", boostBetaZ); // 0.435 for pPb


  // read in from conf file
  int vetoMode    = p.settings.mode("POWHEG:veto");
  int MPIvetoMode = p.settings.mode("POWHEG:MPIveto");
  loadhooks = (vetoMode > 0 || MPIvetoMode > 0);

  if (loadhooks) { // if NOT use SCALUP as starting scale
    if (vetoMode > 0) { // use kinematical limit as starting scale and veto
      p.readString("SpaceShower:pTmaxMatch = 2");
      p.readString("TimeShower:pTmaxMatch = 2");
    }
    if (MPIvetoMode > 0) {
      p.readString("MultipartonInteractions:pTmaxMatch = 2");
    }
    // activate POWHEG compliance
    powhegHooks = new QEDQCDPowhegHooks();
    p.setUserHooksPtr((UserHooks *) powhegHooks);
  }

  
  // variables to keep track of
  //----------------------------------------------------------------------
  int iPhoton; // index of hardest photon in pythia event

  double ptMax, ptTemp; // for searching leading/hardest photon
  bool AreWeightsHistosBooked = false;
  
  //----------------------------------------------------------------------
  fileName = argv[3];
  printf("Showering events in %s\n",fileName.c_str());
  p.readString("Beams:LHEF = " + fileName);
  p.init();

  // skip pythia errors and break, when showering has reached the end of the LHE file
  //----------------------------------------------------------------------
  int counter = 0;
  while (counter++ < 1000000) {
    if (!p.next()) {
      if (p.info.atEndOfFile()) break;
      continue;
    }
      
    p.event.bst(0., 0., boostBetaZ);
    if(h_nEvents->GetBinContent(2) == 1)
      cout << "energy of beam a = " << p.event[1].e() << endl
	   << "energy of beam b = " << p.event[2].e() << endl;
      
    h_nEvents->Fill(1.);
    h_nEvents->Fill(2.);
    h_nEvents->Fill(3.);
    h_nEvents->Fill(4.);

    // only once: read in weight IDs and book histograms for each weight
    //----------------------------------------------------------------------
    if(!AreWeightsHistosBooked){
      ReadInWeightIDs(p, vec_weightsID, isSudaWeight);
      // book histograms for each weight
      for(long unsigned int i = 0; i < vec_weightsID.size(); i++){
	vec_directphoton_pt.push_back( TH1D(Form("h_directphoton_pt_%s",vec_weightsID.at(i).c_str()), "direct photon pt", nPtBins, PtBinsMin, PtBinsMax));
	vec_directphoton_pt_leading.push_back( TH1D(Form("h_directphoton_pt_leading_%s",vec_weightsID.at(i).c_str()), "leading direct photon pt", nPtBins, PtBinsMin, PtBinsMax));
	vec_directphoton_pt_leading_bornveto00.push_back( TH1D(Form("h_directphoton_pt_leading_bornveto00%s",vec_weightsID.at(i).c_str()), "leading direct photon pt bornveto00", nPtBins, PtBinsMin, PtBinsMax));
	vec_directphoton_pt_leading_bornveto20.push_back( TH1D(Form("h_directphoton_pt_leading_bornveto20%s",vec_weightsID.at(i).c_str()), "leading direct photon pt bornveto20", nPtBins, PtBinsMin, PtBinsMax));
	vec_directphoton_pt_leading_bornveto40.push_back( TH1D(Form("h_directphoton_pt_leading_bornveto40%s",vec_weightsID.at(i).c_str()), "leading direct photon pt bornveto40", nPtBins, PtBinsMin, PtBinsMax));
	  
	vec_isodirectphoton_pt.push_back( TH1D(Form("h_isodirectphoton_pt_%s",vec_weightsID.at(i).c_str()), "direct photon pt", nPtBins, PtBinsMin, PtBinsMax));
	vec_isodirectphoton_pt_leading.push_back( TH1D(Form("h_isodirectphoton_pt_leading_%s",vec_weightsID.at(i).c_str()), "leading direct photon pt", nPtBins, PtBinsMin, PtBinsMax));
	  
	vec_chjet_pt.push_back( TH1D(Form("h_chjet_pt_%s",vec_weightsID.at(i).c_str()), "charged jet pt", nPtBins, PtBinsMin, PtBinsMax));
	vec_chjet_pt_leading.push_back( TH1D(Form("h_chjet_pt_leading_%s",vec_weightsID.at(i).c_str()), "leading charged jet pt", nPtBins, PtBinsMin, PtBinsMax));

	vec_isoCone_track_phi.push_back( TH1D(Form("h_isoCone_track_phi_%s",vec_weightsID.at(i).c_str()), "isoCone_track_phi", 100, -1., 1.));
	vec_isoCone_track_eta.push_back( TH1D(Form("h_isoCone_track_eta_%s",vec_weightsID.at(i).c_str()), "isoCone_track_eta", 100, -1., 1.));

	vec_dPhiJetGamma.push_back( TH1D(Form("h_dPhiJetGamma_%s",vec_weightsID.at(i).c_str()), "#Delta #phi_{J#gamma}", 136, -0.1, 3.3));
	vec_dPhiJetGamma_noDeltaPhiCut.push_back( TH1D(Form("h_dPhiJetGamma_noDeltaPhiCut_%s",vec_weightsID.at(i).c_str()), "#Delta #phi_{J#gamma} no cut", 170, -0.1, 3.3));

	vec_xJetGamma.push_back( TH1D(Form("h_xJetGamma_%s",vec_weightsID.at(i).c_str()), "x_{J#gamma} = p_{T}^{Jet} / p_{T}^{#gamma}", 150, 0., 3.));

	vec_chJetTrackMult.push_back( TH1D(Form("h_chJetTrackMult_%s",vec_weightsID.at(i).c_str()), "charged track multiplicity within jets", 50, 0.5, 50.5));

	vec_xObs_pGoing.push_back( TH1D(Form("h_xObs_pGoing_%s",vec_weightsID.at(i).c_str()), "xObs_pGoing", 1000, 0., 0.10));
	vec_xObs_PbGoing.push_back( TH1D(Form("h_xObs_PbGoing_%s",vec_weightsID.at(i).c_str()), "xObs_PbGoing", 1000, 0., 0.10));
	vec_xBjorken_1.push_back( TH1D(Form("h_xBjorken_1_%s",vec_weightsID.at(i).c_str()), "xBjorken 1", 1000, 0., 0.10));
	vec_xBjorken_2.push_back( TH1D(Form("h_xBjorken_2_%s",vec_weightsID.at(i).c_str()), "xBjorken 2", 1000, 0., 0.10));
	vec_xBjorken_1_PDF.push_back( TH1D(Form("h_xBjorken_1_PDF_%s",vec_weightsID.at(i).c_str()), "xBjorken 1 alternative version x1pdf", 1000, 0., 0.10));
	vec_xBjorken_2_PDF.push_back( TH1D(Form("h_xBjorken_2_PDF_%s",vec_weightsID.at(i).c_str()), "xBjorken 2 alternative version x2pdf", 1000, 0., 0.10));

	vec_isoPt.push_back( TH1D(Form("h_isoPt_%s",vec_weightsID.at(i).c_str()), "Pt summed in iso cone", nPtBins, PtBinsMin, PtBinsMax/5.));
	vec_UEPtDensity.push_back( TH1D(Form("h_UEPtDensity_%s",vec_weightsID.at(i).c_str()), "UEPtDensity", nPtBins, PtBinsMin, PtBinsMax/5.));
	vec_isoPt_corrected.push_back( TH1D(Form("h_isoPt_corrected_%s",vec_weightsID.at(i).c_str()), "Pt summed in iso cone", 2*nPtBins, -PtBinsMax/5., PtBinsMax/5.));
	vec_xSecTriggerGamma.push_back(TH1D(Form("h_xSecTriggerGamma_%s",vec_weightsID.at(i).c_str()),"accumulated cross section of trigger photons", 1, -0.5, 0.5));
		  
	AreWeightsHistosBooked = true;
      }
    }
      
    // if Sudakov reweighting is activated, get corresponding weight for this event
    // and reload vector with regular weights * sudaWeight for this event 
    if (isSudaWeight) sudaWeight = p.info.getWeightsDetailedValue("sudakovwgt");
    if(vec_weights.size() != 0) vec_weights.clear();
    for(long unsigned int i = 0; i < vec_weightsID.size(); i++)
      vec_weights.push_back(p.info.getWeightsDetailedValue(vec_weightsID.at(i)) * sudaWeight);
 
    // The actual event analysis starts here.
    ptMax  = -1.;
    ptTemp = -1.;
    iPhoton = -1;
    vector<PseudoJet> vPseudo;
      
    // search for hardest photon in this event
    //----------------------------------------------------------------------
    for (int i = 5; i < p.event.size(); i++) {
      if (p.event[i].id() == 22 && p.event[i].isFinal() && // final photon
	  p.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
	  TMath::Abs(p.event[i].eta()) < etaTPC-jetRadius){// in maximal TPC-minus-iso-cone-radius acceptance
	  
	// find ptMax
	ptTemp = p.event[i].pT();
	if (ptTemp > ptMax) {
	  ptMax = ptTemp;
	  iPhoton = i; // remember index of hardest photon
	}
      }
    }
      
    if(iPhoton > 0){
      // vary born veto to check if enough/too much is cut away
      pyHelp.FillForEachWeight(vec_directphoton_pt_leading_bornveto00, p.event[iPhoton].pT(), vec_weights);
	
      if(ptMax > p.info.getScalesAttribute("uborns")*2.0) h_nEvents->Fill(3.,-1.);
      else pyHelp.FillForEachWeight(vec_directphoton_pt_leading_bornveto20, p.event[iPhoton].pT(), vec_weights);
	
      if(ptMax > p.info.getScalesAttribute("uborns")*3.0) h_nEvents->Fill(4.,-1.);
      else pyHelp.FillForEachWeight(vec_directphoton_pt_leading_bornveto40, p.event[iPhoton].pT(), vec_weights);

      // use following line to ignore events with extreme weights that can cause ugly fluctuations
      // but make sure the cross section does not decrease significantly
      if(ptMax > p.info.getScalesAttribute("uborns")*2.5){
	h_nEvents->Fill(1.,-1.);
	continue; // jump to next event = veto event if hardest photon is 3 times harder than born scale
      }
    }
      
    // set up pseudojets and set up background density in eta band
    //----------------------------------------------------------------------
    for (int i = 5; i < p.event.size(); i++) {
      if (p.event[i].isFinal() && p.event[i].isCharged()) {
	if (TMath::Abs(p.event[i].eta()) < etaTPC){
	  vPseudo.push_back(PseudoJet(p.event[i].px(),p.event[i].py(),p.event[i].pz(),p.event[i].e()));
	}
      }
    }
    cs = new ClusterSequence(vPseudo, jetDef_miguel);
    vJets = sorted_by_pt(cs->inclusive_jets(2.)); // argument: jetminpt
    // loop over charged jets
    //----------------------------------------------------------------------
    if (vJets.size() != 0) {
      for(unsigned int j = 0; j < vJets.size(); j++){
	if(TMath::Abs(vJets.at(j).eta()) > etaTPC-jetRadius) continue;
	pyHelp.FillForEachWeight(vec_chjet_pt, vJets.at(j).pt(), vec_weights);
	if(j == 0)
	  pyHelp.FillForEachWeight(vec_chjet_pt_leading, vJets.at(j).pt(), vec_weights);
      }
    }
	    
    // loop over all direct photons
    //----------------------------------------------------------------------
    for (int i = 5; i < p.event.size(); i++) {
      bool isPhotonIsolated;
      if (p.event[i].id() == 22 && p.event[i].isFinal() && // final photon
	  p.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
	  TMath::Abs(p.event[i].eta()) < etaTPC-jetRadius){       // in maximal TPC-minus-iso-cone-radius acceptance

	// photon as pseudojet for analysis
	PseudoJet photonJet(p.event[i].px(), p.event[i].py(), p.event[i].pz(), p.event[i].e());
	if(photonJet.pt() < 15.) continue;
	if(photonJet.pt() > 30.) continue;
	// calculate ue pt density for a given photon i
	double UEPtDensity = 0.;
	//	  printf("UEPtDensity(p.event, i) = %f\n",UEPtDensity);
	pyHelp.FillForEachWeight(vec_UEPtDensity, UEPtDensity, vec_weights);
	// check isolation
	isPhotonIsolated = pyHelp.IsPhotonIsolatedPowheg(p.event, i, etaTPC-jetRadius, isoConeRadius, isoPtMax, UEPtDensity, vec_isoCone_track_phi, vec_isoCone_track_eta, vec_isoPt, vec_isoPt_corrected, vec_weights);

	// Fill histograms
	//----------------------------------------------------------------------
	pyHelp.FillForEachWeight(vec_directphoton_pt, p.event[i].pT(), vec_weights);
	if(i==iPhoton) pyHelp.FillForEachWeight(vec_directphoton_pt_leading, p.event[i].pT(), vec_weights);

	if(isPhotonIsolated){
	  pyHelp.FillForEachWeight(vec_isodirectphoton_pt, p.event[i].pT(), vec_weights);
	  if(i==iPhoton) pyHelp.FillForEachWeight(vec_isodirectphoton_pt_leading, p.event[i].pT(), vec_weights);	    
	}

	if(vJets.size() > 0 && isPhotonIsolated)
	  for(unsigned int iJet = 0; iJet < vJets.size(); iJet++){
	    bool isJetSeparated = ( TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))) > TMath::Pi()/2. );
	    if(vJets.at(iJet).pt() < 10.) break; // vJets are sorted by pt, break is ok
	    if(iJet == 0) pyHelp.FillForEachWeight(vec_xSecTriggerGamma, 0., vec_weights);  // if there is at least one jet, count trigger photons, but only once 
	    // gamma-jet correlation	 
	    pyHelp.FillForEachWeight(vec_dPhiJetGamma_noDeltaPhiCut, TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))), vec_weights);
	    if(!isJetSeparated) continue;
	    // gamma-jet correlation
	    pyHelp.FillForEachWeight(vec_dPhiJetGamma, TMath::Abs(photonJet.delta_phi_to(vJets.at(iJet))), vec_weights);
	    // x_Jet-gamma
	    vector<PseudoJet> vec_jetConst = vJets.at(iJet).constituents();
	    pyHelp.FillForEachWeight(vec_xJetGamma, vJets.at(iJet).pt()/photonJet.pt(), vec_weights);
	    // charged particle multiplicity in jets
	    pyHelp.FillForEachWeight(vec_chJetTrackMult, vec_jetConst.size(), vec_weights);
	    // x_obs p-going direction
	    pyHelp.FillForEachWeight(vec_xObs_pGoing, pyHelp.XObs_pGoing(vJets.at(iJet), photonJet, p.info.eB()), vec_weights);
	    // x_obs Pb-going direction
	    pyHelp.FillForEachWeight(vec_xObs_PbGoing, pyHelp.XObs_PbGoing(vJets.at(iJet), photonJet, p.info.eB()), vec_weights);
	    // real Bjorken x
	    pyHelp.FillForEachWeight(vec_xBjorken_1, p.info.x1() , vec_weights);
	    pyHelp.FillForEachWeight(vec_xBjorken_2, p.info.x2() , vec_weights);
	    pyHelp.FillForEachWeight(vec_xBjorken_1_PDF, p.info.x1pdf() , vec_weights);
	    pyHelp.FillForEachWeight(vec_xBjorken_2_PDF, p.info.x2pdf() , vec_weights);
	  }

	// print scales of event
	// printf("---------------------------------------\n");
	// printf("p.info.pTHat()   = %f\n", p.info.pTHat());
	// printf("p.info.QFac()    = %f\n", p.info.QFac());
	// printf("p.info.QRen()    = %f\n", p.info.QRen());
	// printf("p.info.scalup()  = %f\n", p.info.scalup());
	// printf("\n");
	  
	  
      } // if direct photon in acceptance
    } // particle for-loop
    delete cs;
  } // end of while loop; break if next file

  // statistics on event generation
  p.stat();

  // write histograms to file ----------------------------------------
  TFile file(rootFileName, "RECREATE");

  h_nEvents->Write();
  
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_directphoton_pt, file, etaTPC-jetRadius); // NB: always "width" scaling applied
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_directphoton_pt_leading, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_directphoton_pt_leading_bornveto00, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_directphoton_pt_leading_bornveto20, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_directphoton_pt_leading_bornveto40, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isodirectphoton_pt, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isodirectphoton_pt_leading, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_chjet_pt, file, etaTPC-jetRadius);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_chjet_pt_leading, file, etaTPC-jetRadius);
  
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isoCone_track_phi, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isoCone_track_eta, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_dPhiJetGamma_noDeltaPhiCut, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_dPhiJetGamma, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xJetGamma, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_chJetTrackMult, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xObs_pGoing, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xObs_PbGoing, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xBjorken_1, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xBjorken_2, file, 1.);
  // pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xBjorken_1_PDF, file, 1.); // second version to call bjorken x gives same results
  // pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xBjorken_2_PDF, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isoPt, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_UEPtDensity, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_isoPt_corrected, file, 1.);
  pyHelp.Add_Histos_Scale_Write2File_Powheg(vec_xSecTriggerGamma, file, 1.);
  
  file.Close();

  if (powhegHooks) delete powhegHooks;
  return 0;

}

//----------------------------------------------------------------------
void ReadInWeightIDs(Pythia &p, vector<string> &vec_weightsID, bool &isSudaWeight){
  
  // check if the sudakov weight from enhanced radiation is present
  for (map<string,double>::iterator it = p.info.weights_detailed->begin();
       it != p.info.weights_detailed->end(); ++it) {
    if (it->first == "sudakovwgt") {
      isSudaWeight = true;
      printf("Sudakov reweighting of hard process is taken into account.\n");
      continue;
    }
  }

  // // if more weights at the same time are used,
  // // e.g. for scale or pdf variation, you can  access them like this
  // for (map<string,double>::iterator it = p.info.weights_detailed->begin();
  //      it != p.info.weights_detailed->end(); ++it) {
  //   if (it->first.find("scales") != std::string::npos){
  //     vec_weightsID.push_back(it->first);
  //   }
  // }

  // insert central value always at first position for convenience 
  for (map<string,double>::iterator it = p.info.weights_detailed->begin();
       it != p.info.weights_detailed->end(); ++it) {
    if (it->first == "central"){ // NB: these strings follow 'lhrwgt_id' in powheg-input.save
      vec_weightsID.insert(vec_weightsID.begin(), it->first);
    }
  }

  printf("Number of weights = %lu\n", vec_weightsID.size());
  for(long unsigned int i = 0; i < vec_weightsID.size(); i++)
    printf("weight description at position %lu: %s\n", i, vec_weightsID.at(i).c_str());

  return;
}

// calculate UE pt density
//----------------------------------------------------------------------
double GetUEPtDensity(Event &event, int iPhoton){
  const double etaBandArea = (0.87*2*0.4) - (0.4*0.4*TMath::Pi()); // (tpc eta acceptance x isocone width) minus isocone area 
  double isoCone_dR = 999.;
  double sumPt = 0.; // reset sum of energy in cone

  for (int iTrack = 5; iTrack < event.size(); iTrack++) {
    if ( !event[iTrack].isFinal() ) continue;
    if ( !event[iTrack].isVisible() ) continue;
    if ( !event[iTrack].isCharged() ) continue;
    if ( TMath::Abs(event[iTrack].eta()) > 0.87 ) continue;
    if ( TMath::Abs(event[iPhoton].eta()) > 0.47 ) continue;

    isoCone_dR = sqrt( pow(pyHelp.CorrectPhiDelta(event[iTrack].phi(), event[iPhoton].phi()), 2)
		       + pow(event[iTrack].eta() - event[iPhoton].eta(), 2) );

    if(isoCone_dR > 0.4 &&
       TMath::Abs(event[iTrack].phi() - event[iPhoton].phi()) < 0.4){
      // printf("isoCone_dR = %f\n", isoCone_dR);
      // printf("TMath::Abs(event[iTrack].phi() - event[iPhoton].phi()) = %f\n", TMath::Abs(event[iTrack].phi() - event[iPhoton].phi()));
      sumPt += event[iTrack].pT();
    }
  }
  
  return sumPt/etaBandArea;
}
