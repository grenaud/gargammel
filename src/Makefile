
CXX      = g++ #-g
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)

CXXFLAGS +=  -std=c++0x -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/src/  -c
#CXXFLAGS +=  -static -std=c++0x -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/src/  -c
LDFLAGS  += ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/PutProgramInHeader.o ${BAMTOOLS}/build/src/api/libbamtools.a  ${LIBGAB}/libgab.a ${LIBGAB}/gzstream/libgzstream.a -lm -lz


all: fragSim deamSim adptSim fasta2fastas damage_patterns2prof mapDamage2prof

adptSim.o:	adptSim.cpp
	${CXX} ${CXXFLAGS} adptSim.cpp

adptSim:	adptSim.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

fragSim.o:	fragSim.cpp
	${CXX} ${CXXFLAGS} fragSim.cpp

fragSim:	fragSim.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

deamSim.o:	deamSim.cpp
	${CXX} ${CXXFLAGS} deamSim.cpp

deamSim:	deamSim.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

damage_patterns2prof.o:	damage_patterns2prof.cpp
	${CXX} ${CXXFLAGS} damage_patterns2prof.cpp

damage_patterns2prof:	damage_patterns2prof.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mapDamage2prof.o:	mapDamage2prof.cpp
	${CXX} ${CXXFLAGS} mapDamage2prof.cpp

mapDamage2prof:	mapDamage2prof.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

fasta2fastas.o:	fasta2fastas.cpp
	${CXX} ${CXXFLAGS} fasta2fastas.cpp

fasta2fastas:	fasta2fastas.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f fragSim.o fragSim deamSim.o deamSim damage_patterns2prof.o damage_patterns2prof mapDamage2prof.o mapDamage2prof adptSim.o adptSim fasta2fastas.o fasta2fastas

