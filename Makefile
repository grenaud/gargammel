SHELL := /bin/bash

OS := $(shell uname)

export BAMTOOLSLIB = $(realpath ./bamtools/build/src)


all: 	src/fragSim src/deamSim src/adptSim src/fasta2fastas art_src_MountRainier/art_illumina_src/art_illumina.o

src/fragSim: libgab/libgab.a bamtools/lib/libbamtools.so
	make -C src

src/fasta2fastas: libgab/libgab.a bamtools/lib/libbamtools.so
	make -C src

src/deamSim: libgab/libgab.a bamtools/lib/libbamtools.so
	make -C src

src/adptSim: libgab/libgab.a bamtools/lib/libbamtools.so
	make -C src

libgab/libgab.h:
	rm -rf libgab/
	git clone --depth 1 --recursive https://github.com/grenaud/libgab.git

libgab/libgab.a: bamtools/lib/libbamtools.so  libgab/libgab.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone  --recursive https://github.com/pezmaster31/bamtools.git # && cd bamtools/ && git reset --hard 2bd8699 # d24d850de17134fe4e7984b26493c5c0a1844b35

bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi  && make
	cp bamtools/build/src/api/bamtools_api_export.h bamtools/src/api && cd ../.. 

# ART is fetched as an upstream tarball and then patched: see
# patches/art_illumina_gargammel.patch for what the patch changes and why.
# The touch after patching stops make from trying to re-run automake/autoconf
# just because Makefile.am and configure.ac are now newer than what they generate.
ARTPATCH = patches/art_illumina_gargammel.patch

art_src_MountRainier/art_illumina_src/art_illumina.o: $(ARTPATCH) #todo: add wget after rm
	rm -rf art_src_MountRainier/ art_src_MountRainier_Linux/ art_src_MountRainier_MacOS/ artsrcmountrainier20160605linuxtgz.tgz artsrcmountrainier20160605macostgz.tgz
ifeq ($(OS),Darwin)
	wget -O artsrcmountrainier20160605macostgz.tgz https://www.dropbox.com/s/6zjipl74de9akg5/artsrcmountrainier2016.06.05macos.tgz?dl=0
	tar xvfz artsrcmountrainier20160605macostgz.tgz
	patch -p1 -d art_src_MountRainier_MacOS/ < $(ARTPATCH)
	touch art_src_MountRainier_MacOS/aclocal.m4 art_src_MountRainier_MacOS/configure art_src_MountRainier_MacOS/config.h.in art_src_MountRainier_MacOS/Makefile.in
	cd art_src_MountRainier_MacOS/ && ./configure && make && cd ..
	ln -s art_src_MountRainier_MacOS  art_src_MountRainier
else
	wget -O artsrcmountrainier20160605linuxtgz.tgz https://www.dropbox.com/s/wf8441vslu1f1nd/artsrcmountrainier20160605linuxtgz.tgz?dl=0
	tar xvfz artsrcmountrainier20160605linuxtgz.tgz
	patch -p1 -d art_src_MountRainier_Linux/ < $(ARTPATCH)
	touch art_src_MountRainier_Linux/aclocal.m4 art_src_MountRainier_Linux/configure art_src_MountRainier_Linux/config.h.in art_src_MountRainier_Linux/Makefile.in
	cd art_src_MountRainier_Linux/ && ./configure && make && cd ..
	ln -s art_src_MountRainier_Linux  art_src_MountRainier
endif

# Statically linked binaries, for copying to a machine that does not have the
# same shared libraries, or for a cluster where they cannot be installed.  This
# builds everything the ordinary way first, then relinks the six programs in
# src/ and art_illumina against the static libc, libstdc++, libz and libgsl.
# libgab and bamtools are linked from their .a archives either way, so only the
# system libraries change.
#
# It needs the static system libraries to be installed, which is a separate
# package from the headers on most distributions: on Debian/Ubuntu that is
# libc6-dev, zlib1g-dev and libgsl-dev, all of which ship the .a alongside the
# .so.  On macOS Apple does not ship a static libc and the link will fail; use
# the ordinary build there.
static: all
	$(MAKE) -C src static
	rm -f art_src_MountRainier/art_illumina
	cd art_src_MountRainier/ && $(MAKE) art_illumina LDFLAGS="-static"
	@echo ""
	@echo "Static binaries:"
	@for b in src/fragSim src/deamSim src/adptSim src/fasta2fastas src/damage_patterns2prof src/mapDamage2prof art_src_MountRainier/art_illumina; do \
		if file $$b 2>/dev/null | grep -q 'statically linked'; then echo "  ok      $$b"; else echo "  DYNAMIC $$b"; fi; \
	done

bacterialex:
	mkdir -p bactDBexample
	cd bactDBexample/ && wget -O clovis.tar.gz  https://www.dropbox.com/s/obmr48d72ahjvhp/clovis.tar.gz?dl=1 && tar xvfz clovis.tar.gz && rm -f clovis.tar.gz  && cd ../
	cd bactDBexample/ && wget -O k14.tar.gz https://www.dropbox.com/s/1pdbqbguw0jfzib/k14.tar.gz?dl=1  && tar xvfz k14.tar.gz && rm -f k14.tar.gz && cd ../

clean:
	make -C libgab clean
	make -C src clean

# Runs each subcomponent on small synthetic inputs and checks its output.
# TESTARGS is handed to the test script, e.g.
#	make test TESTARGS="--only fragSim"
#	make test TESTARGS=--list
test: all
	bash tests/run_tests.sh $(TESTARGS)


.PHONY: all static test
