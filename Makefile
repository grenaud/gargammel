SHELL := /bin/bash



all: 	src/fragSim src/deamSim src/adptSim art_src_MountRainier_Linux/art_illumina_src/art_illumina.o

src/fragSim: libgab/utils.o bamtools/lib/libbamtools.so
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
	git clone --depth 1 --recursive https://github.com/pezmaster31/bamtools.git

bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

art_src_MountRainier_Linux/art_illumina_src/art_illumina.o: #todo: add wget after rm 
	rm -rf art_src_MountRainier_Linux/ artsrcmountrainier20160605linuxtgz.tgz
	wget http://www.niehs.nih.gov/research/resources/assets/docs/artsrcmountrainier20160605linuxtgz.tgz
	tar xvfz artsrcmountrainier20160605linuxtgz.tgz
	cd art_src_MountRainier_Linux/ && ./configure && make && cd ..

clean:
	make -C libgab clean
	make -C src clean


.PHONY: all
