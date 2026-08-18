
CXX      = g++ #-g
BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)

CXXFLAGS +=  -std=c++0x -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/src/  -c
LDFLAGS  += ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/PutProgramInHeader.o ${BAMTOOLS}/build/src/libbamtools.a  ${LIBGAB}/libgab.a ${LIBGAB}/gzstream/libgzstream.a -lm -lz

# STATIC=1 links the binaries against the static libc/libstdc++/libz instead of
# the shared ones, so they can be copied to a machine that has neither.  Only
# the link needs it, the objects are the same either way, which is why this is
# in LDFLAGS and not in CXXFLAGS.  Build it with "make static" from the top
# level rather than setting this by hand, so that art_illumina is built the
# same way.
ifdef STATIC
LDFLAGS  += -static
endif

BINARIES = fragSim deamSim adptSim fasta2fastas damage_patterns2prof mapDamage2prof


all: $(BINARIES)

# The objects are reused, but the binaries have to be relinked, and they are
# already up to date with respect to them after an ordinary build.
static:
	rm -f $(BINARIES)
	$(MAKE) all STATIC=1

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

.PHONY: all static clean

