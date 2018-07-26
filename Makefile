# Include the configuration.
-include Makefile.inc

PYTHIAFLAGS=$(CXX_COMMON) -I$(PREFIX_INCLUDE) -L$(PREFIX_LIB) -lpythia8

ROOT=$(shell /gluster2/h_popp01/software/root6/bin/root-config --cflags --libs)

LHAPDF6=-I$(LHAPDF6_INCLUDE) $(LHAPDF6_LIB)/libLHAPDF.so

MERGE=macro/haddav.C

PYTHIA=src/pythia.cpp

# PYTHIA standalone
pythia:	$(PYTHIA) hendrikshelper.o 
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) $(LHAPDF6) -ldl $(ROOT)

# merge programs
haddav: $(MERGE)
	$(CXX) -o $@ $+ -ldl $(ROOT)

haddav_weightCut: $(MERGE2)
	$(CXX) -o $@ $+ $(PYTHIAFLAGS) -ldl $(ROOT) -I$(INC) $(FASTJET)

# helpful functions for pythia
hendrikshelper.o: src/hendrikshelper.cxx src/hendrikshelper.h
	$(CXX) -c src/hendrikshelper.cxx $(PYTHIAFLAGS) -ldl $(ROOT)
