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
##PYHELPER=../src/PythiaAnalysisHelper.cxx

ShowerAnalysis: $(SHOWER) PythiaAnalysisHelper.o
	$(CXX) $(CXX_COMMON) -o $@ $+ -I../src/ $(PYTHIA8) $(FASTJET) -ldl $(ROOT)

# helpful functions for pythia
PythiaAnalysisHelper.o: ../src/PythiaAnalysisHelper.cxx ../src/PythiaAnalysisHelper.h
	$(CXX) -c ../src/PythiaAnalysisHelper.cxx $(PYTHIA8) $(FASTJET) -ldl $(ROOT)

.PHONY : clean
clean:
	rm -f ShowerAnalysis
	rm -f PythiaAnalysisHelper.o
