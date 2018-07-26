#ifndef _PYTHIAANALYSIS_h_included_
#define _PYTHIAANALYSIS_h_included_

#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "hendrikshelper.h"

Pythia8::Pythia p;

HendriksHelper pyHelp;

int pTHatStartBin = 0;
bool applyPhotonIso = false;

char rootFileName[1024];

bool applyBoost = false;
double boostBetaZ = 0.;

bool MB_veto = false;  

// pthat bin definition from ALICE JJ production at 8 TeV
const int pTHatBins = 18;
double pTHatBin[pTHatBins+1] = { 9.  , 12. , 16. , 21. , 28.,
				 36. , 45. , 57. , 70. , 85.,
				 99. , 115., 132., 150., 169.,
				 190., 212., 235 , 10000. }; 

const double ptMin = 0., ptMax = 300.;
const int ptBins = 300;
const double etaLarge = 3.,
  etaTPC = 0.9,
  etaEMCal = 0.66,
  etaPHOS = 0.12;



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


#endif
