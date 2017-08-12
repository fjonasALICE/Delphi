# Include the configuration.
-include Makefile.inc

PYTHIA=$(CXX_COMMON) -I$(PREFIX_INCLUDE) -L$(PREFIX_LIB) -lpythia8

ROOT=$(shell root-config --cflags --libs)

LHAPDF6=-I$(LHAPDF6_INCLUDE) $(LHAPDF6_LIB)/libLHAPDF.so

MERGE=haddav.C

PYTHIAPHOTON=pythia_photons.cpp

# PYTHIA standalone
pythia_photons: $(PYTHIAPHOTON) hendrikshelper.o 
	$(CXX) -o $@ $+ $(PYTHIA) $(LHAPDF6) -ldl $(ROOT)

# POHWEH showered
powheg_direct_photons: $(POWHEG)
	$(CXX) -o $@  $+  $(PYTHIA) -ldl $(ROOT)

# merge programs
haddav: $(MERGE)
	$(CXX) -o $@ $+ -ldl $(ROOT)

haddav_weightCut: $(MERGE2)
	$(CXX) -o $@ $+ $(PYTHIA) -ldl $(ROOT) -I$(INC) $(FASTJET)

# helpful functions for pythia
hendrikshelper.o: hendrikshelper.cxx hendrikshelper.h
	$(CXX) -c hendrikshelper.cxx $(PYTHIA) -ldl $(ROOT)
