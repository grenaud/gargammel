SHELL := /bin/bash

all: inputfolder/endo/endo.2.fa.fai inputfolder/cont/cont.2.fa.fai

inputfolder:
	mkdir -p inputfolder/
	mkdir -p inputfolder/bact
	mkdir -p inputfolder/endo
	mkdir -p inputfolder/cont


AltaiNea.hg19_1000g.21.mapp.vcf.gz.tbi: inputfolder
	wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
	bedtools intersect -header -a AltaiNea.hg19_1000g.21.mod.vcf.gz  -b 21.bed.gz |bgzip -c > inputfolder/endo/AltaiNea.hg19_1000g.21.mapp.vcf.gz 
	tabix -p vcf inputfolder/endo/AltaiNea.hg19_1000g.21.mapp.vcf.gz
	rm -fv AltaiNea.hg19_1000g.21.mod.vcf.gz

21.fa.gz.fai:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz
	gzip -c -d chr21.fa.gz |sed "s/^>chr21/>21/g" |bgzip -c > 21.fa.gz
	samtools faidx 21.fa.gz

HGDP00521.hg19_1000g.21.mapp.vcf.gz.tbi: inputfolder
	wget http://cdna.eva.mpg.de/denisova/VCF/human/HGDP00521.hg19_1000g.21.mod.vcf.gz
	bedtools intersect -header -a HGDP00521.hg19_1000g.21.mod.vcf.gz  -b 21.bed.gz |bgzip -c > inputfolder/cont/HGDP00521.hg19_1000g.21.mapp.vcf.gz 
	tabix -p vcf inputfolder/cont/HGDP00521.hg19_1000g.21.mapp.vcf.gz 
	rm -fv HGDP00521.hg19_1000g.21.mod.vcf.gz


inputfolder/endo/chr21endo1.fa.gz: 21.fa.gz.fai AltaiNea.hg19_1000g.21.mapp.vcf.gz.tbi
	bcftools consensus -i -f  21.fa.gz  inputfolder/endo/AltaiNea.hg19_1000g.21.mapp.vcf.gz | gzip > inputfolder/endo/chr21endo.gz
	../src/fasta2fastas inputfolder/endo/chr21endo.gz inputfolder/endo/chr21endo
	rm -fv inputfolder/endo/chr21endo.gz
	rm -fv inputfolder/endo/AltaiNea.hg19_1000g.21.mapp.vcf.gz
	rm -fv inputfolder/endo/AltaiNea.hg19_1000g.21.mapp.vcf.gz.tbi

inputfolder/cont/chr21cont1.fa.gz: 21.fa.gz.fai HGDP00521.hg19_1000g.21.mapp.vcf.gz.tbi
	bcftools consensus -i -f  21.fa.gz  inputfolder/cont/HGDP00521.hg19_1000g.21.mapp.vcf.gz  | gzip > inputfolder/cont/chr21cont.gz
	../src/fasta2fastas inputfolder/cont/chr21cont.gz inputfolder/cont/chr21cont
	rm -fv inputfolder/cont/chr21cont.gz
	rm -fv inputfolder/cont/HGDP00521.hg19_1000g.21.mapp.vcf.gz
	rm -fv inputfolder/cont/HGDP00521.hg19_1000g.21.mapp.vcf.gz.tbi


inputfolder/endo/endo.2.fa.fai: inputfolder/endo/chr21endo1.fa.gz
	gzip -c -d  inputfolder/endo/chr21endo1.fa.gz >  inputfolder/endo/endo.1.fa
	gzip -c -d  inputfolder/endo/chr21endo2.fa.gz >  inputfolder/endo/endo.2.fa
	samtools faidx inputfolder/endo/endo.1.fa
	samtools faidx inputfolder/endo/endo.2.fa
	rm -fv inputfolder/endo/chr21endo1.fa.gz
	rm -fv inputfolder/endo/chr21endo1.fa.gz.fai
	rm -fv inputfolder/endo/chr21endo2.fa.gz
	rm -fv inputfolder/endo/chr21endo2.fa.gz.fai


inputfolder/cont/cont.2.fa.fai: inputfolder/cont/chr21cont1.fa.gz
	gzip -c -d  inputfolder/cont/chr21cont1.fa.gz >  inputfolder/cont/cont.1.fa
	gzip -c -d  inputfolder/cont/chr21cont2.fa.gz >  inputfolder/cont/cont.2.fa
	samtools faidx inputfolder/cont/cont.1.fa
	samtools faidx inputfolder/cont/cont.2.fa
	rm -fv inputfolder/cont/chr21cont1.fa.gz
	rm -fv inputfolder/cont/chr21cont1.fa.gz.fai
	rm -fv inputfolder/cont/chr21cont2.fa.gz
	rm -fv inputfolder/cont/chr21cont2.fa.gz.fai

