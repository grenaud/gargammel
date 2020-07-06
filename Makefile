SHELL := /bin/bash


OS := $(shell uname)

all: 	src/fragSim src/deamSim src/adptSim src/fasta2fastas art_src_MountRainier/art_illumina_src/art_illumina.o

src/fragSim: libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

src/fasta2fastas: libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

src/deamSim: libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

src/adptSim: libgab/utils.o bamtools/lib/libbamtools.so
	make -C src

libgab/utils.h:
	rm -rf libgab/
	git clone --depth 1 --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: bamtools/lib/libbamtools.so  libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone  --recursive https://github.com/pezmaster31/bamtools.git && cd bamtools/ && git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi  && make && cd ../..

art_src_MountRainier/art_illumina_src/art_illumina.o: #todo: add wget after rm 
	rm -rf art_src_MountRainier/ art_src_MountRainier_Linux/ art_src_MountRainier_MacOS/ artsrcmountrainier20160605linuxtgz.tgz artsrcmountrainier20160605macostgz.tgz
ifeq ($(OS),Darwin)
	wget -O artsrcmountrainier20160605macostgz.tgz https://www.dropbox.com/s/6zjipl74de9akg5/artsrcmountrainier2016.06.05macos.tgz?dl=0 
	tar xvfz artsrcmountrainier20160605macostgz.tgz
	cd art_src_MountRainier_MacOS/ && ./configure && make && cd ..
	ln -s art_src_MountRainier_MacOS  art_src_MountRainier
else
	wget -O artsrcmountrainier20160605linuxtgz.tgz https://www.dropbox.com/s/wf8441vslu1f1nd/artsrcmountrainier20160605linuxtgz.tgz?dl=0
	tar xvfz artsrcmountrainier20160605linuxtgz.tgz
	cd art_src_MountRainier_Linux/ && ./configure && make && cd ..
	ln -s art_src_MountRainier_Linux  art_src_MountRainier
endif

bacterialex:
	mkdir -p bactDBexample
	cd bactDBexample/ && wget -O clovis.tar.gz  https://www.dropbox.com/s/obmr48d72ahjvhp/clovis.tar.gz?dl=1 && tar xvfz clovis.tar.gz && rm -f clovis.tar.gz  && cd ../
	cd bactDBexample/ && wget -O k14.tar.gz https://www.dropbox.com/s/1pdbqbguw0jfzib/k14.tar.gz?dl=1  && tar xvfz k14.tar.gz && rm -f k14.tar.gz && cd ../

clean:
	make -C libgab clean
	make -C src clean


.PHONY: all
