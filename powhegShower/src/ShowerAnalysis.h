#ifndef _SHOWERANALYSIS_h_included_
#define _SHOWERANALYSIS_h_included_

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
#include "PythiaAnalysisHelper.h"

#include "fastjet/ClusterSequence.hh"

#include "QEDQCDPowhegHooks.h"

using std::cout;
using namespace Pythia8;
using fastjet::PseudoJet;
using fastjet::JetDefinition;
using fastjet::ClusterSequence;
using fastjet::antikt_algorithm;

PythiaAnalysisHelper pyHelp;

const double etaDetector = 0.67;

// jet & iso stuff
const double isoConeRadius = 0.4;
const double isoPtMax=2.0;
const double jetRadius = 0.4;
const JetDefinition jetDef_miguel(antikt_algorithm, jetRadius);


#endif
