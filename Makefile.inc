# Install directory prefixes.
PREFIX_BIN=/gluster2/h_popp01/software/pythia8/bin
PREFIX_INCLUDE=/gluster2/h_popp01/software/pythia8/include
PREFIX_LIB=/gluster2/h_popp01/software/pythia8/lib
PREFIX_SHARE=/gluster2/h_popp01/software/pythia8/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

# FASTJET3 configuration.
FASTJET3_USE=true
FASTJET3_BIN=/gluster2/h_popp01/software/fastjet3/bin/
FASTJET3_INCLUDE=/gluster2/h_popp01/software/fastjet3/include
FASTJET3_LIB=/gluster2/h_popp01/software/fastjet3/lib

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/gluster2/h_popp01/software/hepmc2/bin/
HEPMC2_INCLUDE=/gluster2/h_popp01/software/hepmc2/include
HEPMC2_LIB=/gluster2/h_popp01/software/hepmc2/lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=
LHAPDF5_PLUGIN=LHAPDF5.h

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=/gluster2/h_popp01/software/lhapdf6/bin/
LHAPDF6_INCLUDE=/gluster2/h_popp01/software/lhapdf6/include
LHAPDF6_LIB=/gluster2/h_popp01/software/lhapdf6/lib
LHAPDF6_PLUGIN=LHAPDF5.h

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=
PROMC_LIB=

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=/gluster2/h_popp01/software/root6/bin/
ROOT_INCLUDE=/gluster2/h_popp01/software/root6/include
ROOT_LIB=/gluster2/h_popp01/software/root6/lib

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=
BOOST_LIB=
