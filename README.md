=====================================================================================
  gargammel: simulations of ancient DNA datasets
=====================================================================================

gargammel is a set of programs aimed at simulating ancient DNA fragments. For ancient hominin samples
our program can also simulate various levels of present-day human contamination. It can also add
various levels of bacterial contamination.



Questions :
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
* cmake
* zlib

If you plan on using ms2chromosomes.py to simulate chromosomes based on ms, you also need: 
 * Hudson's ms
 * seq-gen


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
 * endo
 * cont
 * bact

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

  TODO


Example of bacterial databases:
-------------------------------------------------------------------------------------

TODO
