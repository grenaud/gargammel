=====================================================================================
  gargammel: simulations of ancient DNA datasets
=====================================================================================

gargammel is a set of programs aimed at simulating ancient DNA fragments. For ancient hominin samples
our program can also simulate various levels of present-day human contamination. It can also add
various levels of bacterial contamination.



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

Which represent the endogenous ancient human, the present-day human contaminant and the bacterial contamination respectively. Each file inside represents a genome (nor a chromosome or scaffold). The endogenous ancient human can only contain more than 2 genomes since it is a diploid individual. For the bacterial contamination, please add a representative set of bacteria for your sample. Please see the section about the examples of bacterial databases.



Example of usage:
-------------------------------------------------------------------------------------

This is an example of usage to simulate a slightly contaminated (10%) dataset. First we will simulate chromosomes using ms and seq-gen:

    mkdir data
  
Next we will create 1000 simulations of 2 lineages that are allowed to coalesce after 0.2 units of coalescence. The first one will represent our endogenous ancient human while the other, the present-day human contaminant. It will also generate an additiona chromosome from the same population as the contaminant to be used as reference for alignment . We generate sequences for those using the following script:

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

The segsites files corresponds to heterozygous sites between both endogenous genomes.


Then we will create the aDNA fragments:
    cd ..
    ./gargammel.pl -c 3  --comp 0,0.08,0.92 -s src/sizedist.size.gz  -matfile src/matrices/single-  -o data/simulation data/

This will simulated a dataset with 8% human contamination. The deamination rate that will be used will follow a single-strand deamination using the empirical rates measured from the Loschbour individual from:

    Lazaridis, Iosif, et al. "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513.7518 (2014): 409-413.


The size distribution of the aDNA fragments is from a subset of:

    Fu, Qiaomei, et al. "Genome sequence of a 45,000-year-old modern human from western Siberia." Nature 514.7523 (2014): 445-449. 

The read size will be 2x75bp and the Illumina platform being simulated is the HiSeq 2500. The final reads will be found:

    data/out_s1.fq.gz
    data/out_s2.fq.gz


Specifying damage/deamination:
-------------------------------------------------------------------------------------

If you use gargammel.pl or deamSim, you can speficiy deamination/damage using either:

1) Use Briggs model parametes (see Briggs, Adrian W., et al. "Patterns of damage in genomic DNA sequences from a Neandertal." Proceedings of the National Academy of Sciences 104.37 (2007): 14616-14621.)
2) Specify a matrix of deamination rates, we use the following format, the first line is the header:

    	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
    pos	rate_{A->C}	rate_{A->G}	rate_{A->T}	rate_{C->A}	rate_{C->G}	rate_{C->T}	rate_{G->A}	rate_{G->C}	rate_{G->T}	rate_{T->A}	rate_{T->C}	rate_{T->G}

The pos. is the position 0,1... after the fragment beginning/end. The rate is specified using the following: estimate  [estimate_low estimate_high]. For example, 0.3 [0.2 0.4] means that the rate of deamination is 0.3 or 30%.

example of a format:

	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
	0	1.853e-3 [1.726e-3..1.989e-3]	4.064e-3 [3.875e-3..4.263e-3]	3.269e-3 [3.099e-3..3.448e-3]	6.661e-3 [6.254e-3..7.094e-3] 3.057e-3 [2.785e-3..3.355e-3] 8.004e-2 [7.865e-2..8.145e-2] 1.236e-2 [    1.183e-2..1.292e-2] 4.131e-3 [3.828e-3..4.459e-3] 6.703e-3 [6.314e-3..7.116e-3] 3.845e-3 [3.624e-3..4.079e-3] 4.581e-3 [4.339e-3..4.836e-3] 2.169e-3 [2.005e-3..2.347e-3]
	1	1.986e-3 [1.849e-3..2.134e-3]	4.273e-3 [4.070e-3..4.487e-3]	3.030e-3 [2.859e-3..3.211e-3]	5.357e-3 [5.001e-3..5.738e-3] 3.188e-3 [2.916e-3..3.485e-3] 1.427e-2 [1.369e-2..1.488e-2] 9.514e-3 [    9.075e-3..9.974e-3]	3.316e-3 [3.061e-3..3.593e-3] 5.061e-3 [4.743e-3..5.400e-3] 3.421e-3 [3.216e-3..3.639e-3] 4.865e-3 [4.620e-3..5.124e-3]	2.201e-3 [2.038e-3..2.377e-3]

This follows the output of https://bitbucket.org/ustenzel/damage-patterns.git

3) You can use one of the precalculated rates of deamination in src/matrices/. There is a damage from single-strand and a double-strand libraries from the following study: 

    Lazaridis, Iosif, et al. "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513.7518 (2014): 409-413.


Bacterial databases:
-------------------------------------------------------------------------------------

If you simply want to use a uniform probability and do not wish to use a weighted list, if for instance your data is in input/bact/ in fasta files ending with .fa, simply type:

    total=`ls -1  input/bact/*fa |wc -l ` && ls -1 input/bact/*fa  | awk -v total="$total" ' {print $1"\t"(1/total)}' > input/bact/list


Example bacterial databases:
-------------------------------------------------------------------------------------

gargammel comes with these examples of environmental bacterial contamination:

environment1_species:

| genome file               | composition (weight) |
| ------------- | -----:|
| GCF_000421765.1_ASM42176v1 | 0.364962381 |
| GCF_000013745.1_ASM1374v1 | 0.2294623609 |
| GCF_000504265.1_ASM50426v1 | 0.0737186869 |
| GCF_000007625.1_ASM762v1 | 0.0616934061 |
| GCF_000348945.1 | 0.0533369163 |
| GCF_000350565.1 | 0.0409137301 |
| GCF_000332075.2_ASM33207v2 | 0.0360394358 |
| GCF_000010105.1_ASM1010v1 | 0.0301203783 |
| GCF_000257545.3_ASM25754v3 | 0.0275859886 |
| GCF_000146045.2_R64 | 0.0169955819 |
| GCF_000403135.1_ArSv1 | 0.0157299079 |
| GCF_000174915.1_ASM17491v1 | 0.0147518742 |
| GCF_000317145.1_ASM31714v1 | 0.0129647217 |
| GCF_000147015.1_ASM14701v1 | 0.0120021821 |
| GCF_000465325.1 | 0.0097224481 |


environment2_species:

| genome file               | composition (weight) |
| ------------- | -----:|
| GCF_000421145.1_ASM42114v1 | 0.2797168192 |
| GCF_000013605.1_ASM1360v1 | 0.1220485685 |
| GCF_000012145.1_ASM1214v1 | 0.1003648156 |
| GCF_000385945.1_ASM38594v1 | 0.063259174 |
| GCF_000497425.1 | 0.0613174429 |
| GCF_000010105.1_ASM1010v1 | 0.0553500967 |
| GCF_000013745.1_ASM1374v1 | 0.0513526183 |
| GCF_000376425.1_ASM37642v1 | 0.0422420349 |
| GCF_000389985.1_1.0 | 0.0411819471 |
| GCF_000168035.1_ASM16803v1 | 0.0302537136 |
| GCF_000444095.1_ASM44409v1 | 0.0228518329 |
| GCF_000282235.1 | 0.0220871445 |
| GCF_000156475.1_ASM15647v1 | 0.0211813315 |
| GCF_000341815.1 | 0.0170649146 |
| GCF_000007625.1_ASM762v1 | 0.0156512439 |
| GCF_000468475.2_Amad2 | 0.0111628046 |
| GCF_000504265.1_ASM50426v1 | 0.0102003507 |
| GCF_000010605.1_ASM1060v1 | 0.0090874423 |
| GCF_000257545.3_ASM25754v3 | 0.0081422525 |
| GCF_000007025.1_ASM702v1 | 0.0079966465 |
| GCF_000226395.1_PenChr_Nov2007 | 0.0074868053 |
