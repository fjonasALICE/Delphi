# Include the configuration.
-include Makefile.inc

PYTHIAFLAGS=-I$(PYTHIA8)/include -L$(PYTHIA8)/lib/ -L$(PYTHIA8)/lib/archive -lpythia8

ROOT=$(shell root-config --cflags --libs)

LHAPDF6=$(shell lhapdf-config --cflags --ldflags)

FASTJET=$(shell fastjet-config --cxxflags --libs)

MERGE=macro/haddav.C

PYTHIA=src/PythiaAnalysis.cpp
PYTHIATEST=src/PythiaAnalysis.cpp
HELPER=PythiaAnalysisHelper

# PYTHIA standalone
PythiaAnalysis:	$(PYTHIA) $(HELPER).o 
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) $(FASTJET) $(LHAPDF6) -ldl $(ROOT)

# merge programs
haddav: $(MERGE)
	$(CXX) -o $@ $+ -ldl $(ROOT)

haddav_weightCut: $(MERGE2)
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) -ldl $(ROOT) -I$(INC) $(FASTJET)

# helpful functions for pythia
PythiaAnalysisHelper.o: src/PythiaAnalysisHelper.cxx src/PythiaAnalysisHelper.h
	$(CXX) $(CXXCOMMON) -c src/PythiaAnalysisHelper.cxx $(ROOT) $(PYTHIAFLAGS) $(FASTJET) 
