# Include the configuration.
-include Makefile.inc

PYTHIAFLAGS=$(CXX_COMMON) -I$(PREFIX_INCLUDE) -L$(PREFIX_LIB) -lpythia8

ROOT=$(shell /gluster2/h_popp01/software/root6/bin/root-config --cflags --libs)

LHAPDF6=-I$(LHAPDF6_INCLUDE) $(LHAPDF6_LIB)/libLHAPDF.so

MERGE=macro/haddav.C

PYTHIA=src/PythiaAnalysis.cpp
PYTHIATEST=src/PythiaAnalysis.cpp
HELPER=PythiaAnalysisHelper

# PYTHIA standalone
PythiaAnalysis:	$(PYTHIA) $(HELPER).o 
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) $(LHAPDF6) -ldl $(ROOT)

# merge programs
haddav: $(MERGE)
	$(CXX) -o $@ $+ -ldl $(ROOT)

haddav_weightCut: $(MERGE2)
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) -ldl $(ROOT) -I$(INC) $(FASTJET)

# helpful functions for pythia
PythiaAnalysisHelper.o: src/PythiaAnalysisHelper.cxx src/PythiaAnalysisHelper.h
	$(CXX) -c src/PythiaAnalysisHelper.cxx $(PYTHIAFLAGS) -ldl $(ROOT)
