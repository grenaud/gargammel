
CXX      = g++ #-g
LIBGAB   = /home/gabriel/install/lib/libgab/

CXXFLAGS =  -std=c++0x -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/  -c
LDFLAGS  = -lz


all: extractFragments 

extractFragments.o:	extractFragments.cpp
	${CXX} ${CXXFLAGS} extractFragments.cpp


extractFragments:	extractFragments.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f extractFragments.o extractFragments

