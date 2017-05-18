# Include the configuration.
-include Makefile.inc

PYTHIA=-I$(PREFIX_INCLUDE) -L$(PREFIX_LIB) -lpythia8

ROOT=$(shell root-config --cflags --libs)

MERGE=haddav.C

PYTHIAPHOTON=pythia_photons.cpp

# PYTHIA standalone
pythia_photons: $(PYTHIAPHOTON) 
	$(CXX) -o $@ $+ $(PYTHIA) -ldl $(ROOT)


# POHWEH showered
powheg_direct_photons: $(POWHEG)
	$(CXX) -o $@  $+  $(PYTHIA) -ldl $(ROOT)

# merge programs
haddav: $(MERGE)
	$(CXX) -o $@ $+ -ldl $(ROOT)

haddav_weightCut: $(MERGE2)
	$(CXX) -o $@ $+ $(PYTHIA) -ldl $(ROOT) -I$(INC) $(FASTJET)

