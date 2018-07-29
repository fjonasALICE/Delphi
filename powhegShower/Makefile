ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC

FASTJET=$(shell /lustre/nyx/alice/users/hpoppenb/fastjet-3.3.1-install/bin/fastjet-config --cxxflags --libs)

PYTHIADIR=/lustre/nyx/alice/users/hpoppenb/pythia8226-install
PYTHIAINC=$(PYTHIADIR)/include
PYTHIALIB=$(PYTHIADIR)/lib
PYTHIA8=-I$(PYTHIAINC) -L$(PYTHIALIB) -lpythia8

ROOT=$(shell root-config --cflags --libs)

SHOWER=src/ShowerAnalysis.cpp
HELPER=ShowerAnalysisHelper

ShowerAnalysis: $(SHOWER)
	$(CXX) $(CXX_COMMON) -o $@ $+ $(PYTHIA8) $(FASTJET) -ldl $(ROOT)

# helpful functions for pythia
ShowerAnalysisHelper.o: ShowerAnalysisHelper.cxx ShowerAnalysisHelper.h
	$(CXX) -c ShowerAnalysisHelper.cxx $(PYTHIA) -ldl $(ROOT)

.PHONY : clean
clean:
	rm -f ShowerAnalysis
	rm -f ShowerAnalysisHelper.o
