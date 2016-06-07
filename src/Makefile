
CXX      = g++ #-g
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)

CXXFLAGS =  -std=c++0x -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -c
LDFLAGS  = ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/utils.o ${LIBGAB}/gzstream/libgzstream.a -lm -lz


all: fragSim deamSim

fragSim.o:	fragSim.cpp
	${CXX} ${CXXFLAGS} fragSim.cpp

fragSim:	fragSim.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

deamSim.o:	deamSim.cpp
	${CXX} ${CXXFLAGS} deamSim.cpp

deamSim:	deamSim.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f fragSim.o fragSim deamSim.o deamSim

