g++ ${1} \
-O2 -ansi -W -Wall -std=c++11 -Wshadow -m64 -Wno-shadow \
-o ${1}.exe \
-I$PYTHIA8/include -L$PYTHIA8/lib/ -L$PYTHIA8/lib/archive -lpythia8 \
-L/home/florianjonas/tools/alice/sw/ubuntu1804_x86-64/GMP/latest/lib -l:libgmp.a \
`root-config --cflags --ldflags --glibs ` \
`fastjet-config --cxxflags --libs` \
`lhapdf-config --cflags --ldflags` \