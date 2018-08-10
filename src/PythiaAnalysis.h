#ifndef _PYTHIAANALYSIS_h_included_
#define _PYTHIAANALYSIS_h_included_

#include <iostream>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "PythiaAnalysisHelper.h"

#include "fastjet/ClusterSequence.hh"

using std::cout;
using namespace Pythia8;
using fastjet::PseudoJet;
using fastjet::JetDefinition;
using fastjet::ClusterSequence;
using fastjet::antikt_algorithm;

Pythia8::Pythia p;

PythiaAnalysisHelper pyHelp;

int pTHatStartBin = 0; // option to skip the first pthat bins
bool producePhotonIsoSpectra = false;
bool useGammaJetCorrelations = false;

char rootFileName[1024]; // output file name

bool applyBoost = false;
double boostBetaZ = 0.; // boost in z direction, needed for asymmetric collision systems

bool MB_veto = false; // true -> omit double counting between pthatbins and MB generation


// kinematic range
//const double ptMin = 0., ptMax = 300.;
//const int ptBins = 97;


const double yDefault = 0.8,
  etaLarge = 3.,
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


// jet & iso stuff
const double isoConeRadius = 0.4;
const double isoPtMax=1.5;
const double jetRadius = 0.4;
const JetDefinition jetDef_miguel(antikt_algorithm, jetRadius);

double photonPtMax;
double photonPtTemp;
int iPhoton;


#endif
