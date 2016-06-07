SHELL := /bin/bash



all: 	src/fragSim

src/fragSim: libgab/utils.o bamtools/lib/libbamtools.so
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

clean:
	make -C libgab clean
	make -C src clean


.PHONY: all
