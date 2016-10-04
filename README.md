=====================================================================================
  gargammel: simulations of ancient DNA datasets
=====================================================================================

gargammel is a set of programs aimed at simulating ancient DNA fragments. For ancient hominin samples
our program can also simulate various levels of present-day human contamination and microbial contamination.

The website for gargammel can be found here: https://grenaud.github.io/gargammel/


Questions/feature requests :
-------------------------------------------------------------------------------------

	contact: Gabriel Renaud   
	email:	 gabriel [dot] reno [ at sign ] gmail.com


Downloading:
-------------------------------------------------------------------------------------
Do a :

    git clone --recursive  --depth 1 https://github.com/grenaud/gargammel.git


Requirements:
-------------------------------------------------------------------------------------
* git
* c++ compiler supporting C11++
* cmake, you can install on Ubuntu by typing: sudo apt install cmake
* zlib
* lib gsl, you can install on Ubuntu by typing: sudo apt-get install libgsl0-dev

If you plan on using ms2chromosomes.py to simulate chromosomes based on ms, you also need: 
 * Hudson's ms (see: http://home.uchicago.edu/rhudson1/source/mksamples.html)
 * seq-gen, you can install on Ubuntu by typing:   sudo apt install seq-gen

Both should be installed in your path. 

Installation:
-------------------------------------------------------------------------------------

In the main directory, simply type

  make 

This should install bamtools (C++ library to read/write BAM files) and ART (Illumina read simulator).


Overview:
-------------------------------------------------------------------------------------

The main driver script, gargammel.pl calls the following programs in order to 
simulate the in vivo process by which ancient DNA fragments are retrieved:

* fragSim: simulation of ancient DNA fragments being retrieved at random from the genome
* deamSim: simulation of damage to the fragments selected by fragSim
* adptSim: adding of adapters to create raw Illumina reads (without errors and quality scores)

Finally, the simulated raw Illumina reads are sent to ART to add sequencing errors and corresponding quality scores.

Input description:
-------------------------------------------------------------------------------------

The basic input is a directory with 3 subfolders named:
 * endo/
 * cont/
 * bact/

Which represent the endogenous ancient human, the present-day human contaminant and the microbial contamination respectively. Each file inside represents a genome (nor a chromosome or scaffold). The endogenous ancient human can only contain more than 2 genomes since it is a diploid individual. For the microbial contamination, please add a representative set of microbes for your sample (see the section about the examples of microbial databases).



Example of usage:
-------------------------------------------------------------------------------------

This is an example of usage to simulate a slightly contaminated (8%) dataset. First, we will simulate chromosomes using ms and seq-gen:

    mkdir data
  
Next, we will create 1000 simulations of 2 lineages that are allowed to coalesce after 0.2 units of coalescence. The first one will represent our endogenous ancient human while the other, the present-day human contaminant. It will also generate an additional chromosome from the same population as the contaminant to be used as reference for alignment. We generate sequences for those using the following script:

    cd data/
    python ../ms2chromosomes.py  -s 0.2 -f . -n 1000 
    rm -rfv simul_* seedms #cleanup
  
This will create the following files:

    cont/cont.0.fa
    cont/cont.1.fa
    endo/endo.1.fa
    endo/endo.2.fa
    endo/segsites
    ref.fa

The segsites files correspond to heterozygous sites between both endogenous genomes.


Then we will create the aDNA fragments:

    cd ..
    ./gargammel.pl -c 3  --comp 0,0.08,0.92 -f src/sizefreq.size.gz  -matfile src/matrices/single-  -o data/simulation data/

This will simulate a dataset with 8% human contamination. The deamination rate that will be used will follow a single-strand deamination using the empirical rates measured from the Loschbour individual from:

    Lazaridis, Iosif, et al. "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513.7518 (2014): 409-413.


The size distribution of the aDNA fragments is a subset of:

    Fu, Qiaomei, et al. "Genome sequence of a 45,000-year-old modern human from western Siberia." Nature 514.7523 (2014): 445-449. 

The read size will be 2x75bp and the Illumina platform being simulated is the HiSeq 2500. The final reads will be found:

    data/out_s1.fq.gz
    data/out_s2.fq.gz


Here are further examples of usage:

* Low coverage 0.5X coverage with fragments of length 40:

`gargammel.pl -c 0.5  --comp 0,0,1 -l 40    -o data/simulation data/`

* Generating exactly 1M fragments of length with a log-normal distribution of location 4.106487474 and scale 0.358874723:

`gargammel.pl -n 1000000  --comp 0,0,1 --loc  4.106487474 --scale  0.358874723   -o data/simulation data/`

* High coverage (20X) with high amount of present-day contamination (40%) with fragments of length 45:

`gargammel.pl -c 20  --comp 0,0.4,0.6 -l 45 -o data/simulation data/`

* Evaluating the impact of mapping 1M fragments with length 40 without double-stranded deamination:

`gargammel.pl -n 1000000  --comp 0,0,1 -l 40    -o data/simulation data/`

* Evaluating the impact of mapping 1M fragments    with length 40 with double-stranded deamination:

`gargammel.pl -n 1000000  --comp 0,0,1 -l 40 -damage 0.03,0.4,0.01,0.3   -o data/simulation data/`

* Generate a single-end run of 96 cycles on a HiSeq 2500 Illumina run with 1M fragments of 40bp:

`gargammel.pl -n 1000000  --comp 0,0,1 -l 40 -rl 96  -se -ss HS25 -o data/simulation data/`

* Generate a paired-end run of 96 cycles on a HiSeq 2500 Illumina run with 1M fragments of 40bp:

`gargammel.pl -n 1000000  --comp 0,0,1 -l 40 -rl 96      -ss HS25 -o data/simulation data/`





Specifying damage/deamination:
-------------------------------------------------------------------------------------

If you use gargammel.pl or deamSim, you can speficiy deamination/damage using either:

1. Use Briggs model parametes (see Briggs, Adrian W., et al. "Patterns of damage in genomic DNA sequences from a Neandertal." Proceedings of the National Academy of Sciences 104.37 (2007): 14616-14621.)
2. Specify a matrix of deamination rates, we use the following format, the first line is the header:

    	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
    	pos	rate_{A->C}	rate_{A->G}	rate_{A->T}	rate_{C->A}	rate_{C->G}	rate_{C->T}	rate_{G->A}	rate_{G->C}	rate_{G->T}	rate_{T->A}	rate_{T->C}	rate_{T->G}

The pos. is the position 0,1... after the fragment beginning/end. The rate is specified using the following: estimate  [estimate_low estimate_high]. For example, 0.3 [0.2 0.4] means that the rate of deamination is 0.3 or 30%.

example of a format:

	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
	0	1.853e-3 [1.726e-3..1.989e-3]	4.064e-3 [3.875e-3..4.263e-3]	3.269e-3 [3.099e-3..3.448e-3]	6.661e-3 [6.254e-3..7.094e-3] 3.057e-3 [2.785e-3..3.355e-3] 8.004e-2 [7.865e-2..8.145e-2] 1.236e-2 [    1.183e-2..1.292e-2] 4.131e-3 [3.828e-3..4.459e-3] 6.703e-3 [6.314e-3..7.116e-3] 3.845e-3 [3.624e-3..4.079e-3] 4.581e-3 [4.339e-3..4.836e-3] 2.169e-3 [2.005e-3..2.347e-3]
	1	1.986e-3 [1.849e-3..2.134e-3]	4.273e-3 [4.070e-3..4.487e-3]	3.030e-3 [2.859e-3..3.211e-3]	5.357e-3 [5.001e-3..5.738e-3] 3.188e-3 [2.916e-3..3.485e-3] 1.427e-2 [1.369e-2..1.488e-2] 9.514e-3 [    9.075e-3..9.974e-3]	3.316e-3 [3.061e-3..3.593e-3] 5.061e-3 [4.743e-3..5.400e-3] 3.421e-3 [3.216e-3..3.639e-3] 4.865e-3 [4.620e-3..5.124e-3]	2.201e-3 [2.038e-3..2.377e-3]

This follows the output of https://bitbucket.org/ustenzel/damage-patterns.git

3. You can use one of the precalculated rates of deamination in src/matrices/. There is a damage from single-strand and a double-strand libraries from the following study: 

    Lazaridis, Iosif, et al. "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513.7518 (2014): 409-413.


How can I get an ancient misincorporation profile for gargammel?
-------------------------------------------------------------------------------------

You could generate it manually, the format is as follows:

    # comment
    Chr	End	Std	Pos	A	C	G	T	Total
    [chr]	['5p' or '3p']	['+' or '-']	[pos wrt the 5p/3p end]	[count A]	[count C]	[count G]	[count T]	[sum of counts]
 
For instance:

	# table produced by mapDamage version 2.0.5-1-ge06bd84
	# using mapped file Ust_Ishim.hg19_1000g.bam and human_g1k_v37.fasta as reference file
	# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads
	Chr	End	Std	Pos	A	C	G	T	Total
	21	3p	+	-4	177086	83624	114115	150943	525768
	21	3p	+	-3	191241	80099	104155	150269	525764
	21	3p	+	-2	197747	63995	127660	136360	525762
	21	3p	+	-1	180637	49770	79519	215833	525759
	21	3p	+	1	188505	79678	204246	53417	525846
	21	3p	+	2	156848	74009	128222	166767	525846
	21	3p	+	3	188608	75382	113613	148243	525846
	21	3p	+	4	173245	84205	117226	151170	525846

The lines above specify the base count close +/- 4 bases to the 3p end for fragments mapping to the + strand. An example of this type of file is found here: src/dnacomp.txt

Such a file can be generated using mapDamage2.0: 

	Jónsson, Hákon, et al. "mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters." Bioinformatics (2013): btt193.

It is normally called "dnacomp.txt" in the output directory, you can filter a single chromosome (in this case 21) using this command:

	grep "^21\|^#\|^Chr"  /path to mapDamage output/results_[sample name]/dnacomp.txt >  dnacomp.txt


How can I specify the size distribution?
-------------------------------------------------------------------------------------

Ancient DNA molecules tend to be fragmented and can be very short but tend to have a specific shape. Both for the wrapper script (gargammel.pl) and the fragment simulation program (fragSim), there are are 4 ways to specifiy the :
1) Specify a fixed length using -l 
2) Open a file containing the size distribution using -s, one empirical fragment length per line eg:
------------
82
95
66
144
87
68
74
48
77
43
------------
3) Open a file containing the size frequencies using -f in the format "size[TAB]freq" eg:
------------
40	0.017096
41	0.01832
42	0.0201954
43	0.018399
44	0.0195637
45	0.0198993
46	0.0196822
47	0.0209456
48	0.0203929
49	0.0199783
50	0.0204323
------------
4) Specify the size distribution using parameters from a log-normal distribution, using options --loc and --scale.

How can I get parameters for the size distribution?
-------------------------------------------------------------------------------------

If you wish to specify the aDNA fragment size distribution as a log-normal, you can use the following script to infer the location and scale parameters:

	#!/usr/bin/env Rscript-3.2.0
	library(fitdistrplus)
	library(MASS)
		
	args=(commandArgs(TRUE))
	
	data <- read.table(args[1]);
		
	df<-fitdistr(data$V1, "lognormal")
	
	print(df);

You can change the header to suit the version of R that you have.

Bacterial databases:
-------------------------------------------------------------------------------------


For the input/bact/ directory which represent the microbial contamination, gargammel needs a set of fasta files that represent the different microbes. Each file corresponds to exactly one microbial species. Each fasta file must contain the genome of the microbial species, multiple scaffolds and plasmids are allowed. Each fasta file must also be faidx indexed. This directory must also contain a file called "list".  This file contains the list of every fasta files in that directory along with their relative abundance in the desired bacterial contamination. For example:

    bacteria1.fa	0.5
    bacteria2.fa	0.3
    bacteria3.fa	0.2


Examples of bacterial databases:
-------------------------------------------------------------------------------------

If you wish to download an example of a suitable bacterial database, you can simply type:
   
     make bacterialex

this will create a directory called bactDBexample/ which contains clovis/ and k14/, the profiled microbial communities from Rasmussen et al. "The genome of a Late Pleistocene human from a Clovis burial site in western Montana." Nature 506.7487 (2014): 225-229. and Seguin-Orlando et al. "Genomic structure in Europeans dating back at least 36,200 years." Science 346.6213 (2014): 1113-1118, respectively.

You can copy the files from the fasta/ directory into the input's bact/ directory as such
    
    cp -v bacterialex/clovis/fasta/* [path to input]/bact/

Creating bacterial databases from a metaBIT:
-------------------------------------------------------------------------------------

metaBIT [https://bitbucket.org/Glouvel/metabit] is a metagenomic profiler from high-throughput sequencing shotgun data. To download the fasta files based on a profile obtained using metaBIT's output, simply supply the "all_taxa.tsv" file, which details the different species and their abundances, make sure you are connected to the internet and use the retrieveFromMetabit script in as such:

    mkdir exampleBacteriaDB
    cd exampleBacteriaDB
    [copy the all_taxa.tsv in the current directory]
    src/microbial_fetcher/retrieveFromMetabit all_taxa.tsv

If you wish, you can enter your email for the ftp from NCBI (to avoid getting banned from the FTP):

   src/microbial_fetcher/retrieveFromMetabit all_taxa.tsv anonymous@server.net


This will download the necessary files from NCBI to create a database suitable for gargammel to simulate microbial species in the exampleBacteriaDB/fasta and run samtools faidx on each file. You need standard UNIX utilities such as awk/sed/python/curl/wget/gzip to be installed as well as samtools. Please move the fasta/ directory produced (exampleBacteriaDB/fasta in the example above) to the input/bact/. The file named "exampleBacteriaDB/fastafasta/list" is the list of bacterial species along with their abundance. Another file, "exampleBacteriaDB/Microbial_ID.log" details the strain/ID and ftp link used.

If you want to use a uniform probability instead of a weighted list, go to "input/bact" and type (if fasta files end with .fa):

    total=`ls -1  input/bact/*fa |wc -l ` && ls -1 input/bact/*fa  | awk -v total="$total" ' {print $1"\t"(1/total)}' > input/bact/list


metaBIT ref: Louvel et al. "metaBIT, an integrative and automated metagenomic pipeline for analyzing microbial profiles from high-throughput sequencing shotgun data." Molecular ecology resources (2016).




Tutorial using empirical sequences for simulations:
-------------------------------------------------------------------------------------

To provide an example of using empirical VCF files to create sequences for the simulation, there is a folder exampleSeq/ with a Makefile. This makefile provides a simple example of creating 2 chromosomes (2 endogenous sequences + 2 contaminant sequences for a diploid genome) from VCF files. This makefile needs the following commands to be installed in the path:

* bedtools
* bgzip
* tabix
* samtools
* bcftools, must support "consensus" command

Make sure that you are connected to the internet and type:

    cd  exampleSeq/ 
    make

This will download the VCF files from the Altai Neanderthal (endogenous) and a present-day human of European descent (contaminant) create 4 files:
         
    inputfolder/endo/endo.2.fa
    inputfolder/endo/endo.1.fa
    inputfolder/cont/cont.1.fa
    inputfolder/cont/cont.2.fa

along with their respective fasta index. If you wish to add bacterial sequences to the mix, please see the section above about "Examples of bacterial databases" and you can copy some files to cont/ directory: cp -v  ../bactDBexample/k14/fasta/* inputfolder/bact/ 

To create a sample with say 10% present-day human contamination with fragment length of 40bp, run:
      
./gargammel.pl -c 0.5  --comp 0,0.1,0.9 -l 40    -o exampleSeq/simulationc10 exampleSeq/inputfolder/

If you have some microbial sequences, to create a sample with say 70% bacterial content, 5% present-day human contamination and 25% endogenous, run:
      
./gargammel.pl -c 0.5  --comp 0.7,0.05,0.25 -l 40    -o exampleSeq/simulationb70c5 exampleSeq/inputfolder/

