
  gargammel: simulations of ancient DNA datasets
=====================================================================================

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gargammel/README.html)


gargammel is a set of programs aimed at simulating ancient DNA fragments. For ancient hominin samples
our program can also simulate various levels of present-day human contamination and microbial contamination.

The website for gargammel can be found here: https://grenaud.github.io/gargammel/


Questions/bug report/feature requests :
-------------------------------------------------------------------------------------

If you have Github account, consider creating an issue, you will help others who might have the same problem.

	contact: Gabriel Renaud   
	email:	 gabriel [dot] reno [ at sign ] gmail.com

I accept pull request for novel features.

Downloading:
-------------------------------------------------------------------------------------
Do a :

    git clone --recursive  --depth 1 https://github.com/grenaud/gargammel.git

or via (bio)conda

```bash
conda install -c bioconda gargammel
```

> Installing with conda will only provide the main gargammel program, for the additional scripts in the repository, please run `git clone` as above, and create the conda environment described below.


Requirements:
-------------------------------------------------------------------------------------
* git
* C++ compiler supporting C++11
* cmake, you can install on Ubuntu by typing: sudo apt install cmake
* zlib
* lib gsl, you can install on Ubuntu by typing: sudo apt-get install libgsl0-dev

If you plan on using ms2chromosomes.py to simulate chromosomes based on ms, you also need: 
 * Hudson's ms (see: http://home.uchicago.edu/rhudson1/source/mksamples.html)
 * seq-gen, you can install on Ubuntu by typing:   sudo apt install seq-gen

Both should be installed in your path.

Alternatively, you can use the supplied [conda](https://https://conda.io/) `environment.yml` file to download and set up all dependencies described in this README for you.


    conda env create -f environment.yml


Installation:
-------------------------------------------------------------------------------------

> If you are using the conda enviroment, you can skip this step and just load the environment with `conda activate gargammel`. All subsequent steps you can replace `gargammel.pl` with just `gargammel`.

In the main directory, simply type

  make 

This should install bamtools (C++ library to read/write BAM files) and ART (Illumina read simulator).

### Static binaries

If you need binaries that can be copied to another machine, or to a cluster
where you cannot install libraries, type

  make static

instead. This builds everything as above, then relinks the six programs in
`src/` and `art_illumina` against the static libc, libstdc++, libz and libgsl,
so that they have no shared library dependencies at all:

    $ ldd src/fragSim
    	not a dynamic executable

The seven binaries are then self-contained; `gargammel.pl` itself is a Perl
script and still needs perl. The simulated reads are unaffected: a static build
gives the same output as an ordinary one, seed for seed.

This needs the static system libraries, which on most distributions are a
separate package from the headers. On Debian/Ubuntu they come with `libc6-dev`,
`zlib1g-dev` and `libgsl-dev`. The linker warns that `getaddrinfo` in a static
binary needs the glibc it was linked against; this comes from a part of
bamtools that gargammel never calls and can be ignored. macOS does not ship a
static libc, so use the ordinary build there.

To go back to ordinary dynamic binaries, run `make clean` first: `make` on its
own will leave the static ones in place, since they are newer than the sources.

Tests:
-------------------------------------------------------------------------------------

To check that the individual subcomponents behave as expected, type

  make test

This runs fragSim, deamSim, adptSim, fasta2fastas, mapDamage2prof and
damage_patterns2prof on small synthetic inputs, then runs gargammel.pl itself,
and verifies the output of each. It takes about a minute, writes everything to
a temporary directory and needs no data beyond what is in this repository.

To run a single group of tests, or to keep the temporary directory around to
look at what was produced:

  make test TESTARGS="--only fragSim"
  bash tests/run_tests.sh --list
  bash tests/run_tests.sh --only 'deamSim|adptSim' --keep

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

Which represent the endogenous ancient human, the present-day human contaminant and the microbial contamination respectively. Each file inside represents a genome (not simply a chromosome or scaffold). The endogenous ancient human can only contain more than 2 genomes since it is a diploid individual. For the microbial contamination, please add a representative set of microbes for your sample (see the section about the examples of microbial databases).



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

This will simulate a dataset with 8% human contamination. The rate of misincorporation due to deamination that will be used will follow a single-strand deamination using the empirical rates measured from the Loschbour individual from:

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

* Reproduce a simulation exactly, by fixing the seed:

`gargammel.pl -n 1000000  --comp 0,0,1 -l 40 -rl 96      -ss HS25 --seed 31337 -o data/simulation data/`



Simulating ancient DNA in one line:
-------------------------------------------------------------------------------------

gargammel.pl exists to mix endogenous, present-day human and microbial sources
in a given proportion. If you only have one genome to draw from, the four
programs it drives can be chained directly: each stage reads what the previous
one writes, so the whole simulation is a single command that never touches the
disk in between:

    src/fragSim -n 1000000 -f src/sizefreq.size.gz --seed 1 ref.fa \
     | src/deamSim -damage 0.024,0.36,0.0097,0.68 --seed 2 /dev/stdin \
     | src/adptSim -l 75 -artp /dev/stdout --seed 3 /dev/stdin \
     | art_src_MountRainier/art_illumina -ss HS25 -amp -na -p -l 75 -c 1 -rs 4 \
        -i - --fq1 sim_s1.fq.gz --fq2 sim_s2.fq.gz

That is the whole pipeline: fragSim draws the fragments, deamSim damages them,
adptSim adds the adapters and splits each fragment into the two mates, and
art_illumina adds the sequencing errors and the quality scores. The result is
`sim_s1.fq.gz` and `sim_s2.fq.gz`, exactly what gargammel.pl would have
produced for a single source.

`ref.fa` needs a `ref.fa.fai` next to it, which `samtools faidx ref.fa`
creates. Give each program its own `--seed` (`-rs` for art) if you want the run
to be reproducible; see the section on reproducible simulations. deamSim and
adptSim read their input as `/dev/stdin` rather than `-`.

For a single-end run, ask adptSim for one read per fragment and drop `-p`:

    src/fragSim -n 1000000 -l 45 --seed 1 ref.fa \
     | src/deamSim -damage 0.024,0.36,0.0097,0.68 --seed 2 /dev/stdin \
     | src/adptSim -l 75 -arts /dev/stdout --seed 3 /dev/stdin \
     | art_src_MountRainier/art_illumina -ss HS25 -amp -na -l 75 -c 1 -rs 4 \
        -i - -o - | gzip -c > sim_s.fq.gz

### Keeping the intermediate files with tee

Piping everything means the intermediate stages are gone by the time the reads
are written, and those are usually what you want to compare the reads against:
the undamaged fragments tell you what each read should have been, and the
deflines carry the true coordinates. Add a `tee` between any two stages to keep
one without interrupting the flow:

    src/fragSim -n 1000000 -f src/sizefreq.size.gz --seed 1 ref.fa \
     | tee >(gzip -c > sim.e.fa.gz) \
     | src/deamSim -damage 0.024,0.36,0.0097,0.68 --seed 2 /dev/stdin \
     | tee >(gzip -c > sim_d.fa.gz) \
     | src/adptSim -l 75 -artp /dev/stdout --seed 3 /dev/stdin \
     | tee >(gzip -c > sim_a.fa.gz) \
     | art_src_MountRainier/art_illumina -ss HS25 -amp -na -p -l 75 -c 1 -rs 4 \
        -i - --fq1 sim_s1.fq.gz --fq2 sim_s2.fq.gz

This writes the same files gargammel.pl leaves behind (`.e.fa.gz` the
fragments, `_d.fa.gz` the damaged fragments, `_a.fa.gz` the amplicons handed to
art) and the reads are unchanged: a `tee` only copies the stream. Keep the ones
you need and drop the rest; comparing `sim.e.fa.gz` with `sim_d.fa.gz` position
by position, for instance, gives you the exact set of deaminated bases.

`>(...)` is a bash process substitution, so run these under bash rather than a
plain POSIX shell. If your shell does not have it, `tee sim_d.fa` and gzipping
afterwards does the same thing at the cost of the uncompressed file. Note that
the shell reports the exit status of the last command in a pipeline only; use
`set -o pipefail` if you want a failure in fragSim or deamSim to be noticed.


Reproducible simulations:
-------------------------------------------------------------------------------------

Pass `--seed [int]` to gargammel.pl and two runs with the same seed and the
same arguments produce byte-identical output. Without it, every run is seeded
from the clock, as before.

The seed is spread over the whole pipeline: gargammel.pl seeds its own draws
(it splits the requested number of fragments over the input genomes at random)
and hands each of fragSim, deamSim, adptSim and art a separate seed derived
from it, so the individual programs do not share a random stream. Each of those
programs also accepts its own `--seed` (`-rs` for art) if you drive them
directly.

All four sub-programs had to be seedable for this to hold. adptSim in
particular pads fragments shorter than the read length with random bases, and
that padding used to come from a clock-seeded generator.

art_illumina now carries its own random number generator rather than using the
one from the C library. It reproduces glibc's `rand()` exactly, so runs match
earlier results on Linux, and art_illumina on its own gives the same reads from
the same seed on Linux and on macOS.

Note that this cross-platform guarantee covers art_illumina only. fragSim,
deamSim and adptSim still draw through the C library (`rand`, `drand48`) and
`std::default_random_engine`, all of which differ between implementations, so a
whole gargammel run with a fixed seed reproduces exactly on the same machine and
toolchain, not between Linux and macOS.


Speed and I/O of the read simulator:
-------------------------------------------------------------------------------------

ART is the slowest stage of the pipeline. The copy gargammel builds is patched
(see `patches/art_illumina_gargammel.patch`, applied by the Makefile when ART is
unpacked) and runs several times faster than the stock version on
gargammel-style input. On 200,000 amplicons drawn with the size distribution in
`src/sizefreq.size.gz`, paired-end at 75bp, it takes 0.9s against 7.3s;
single-end at 75bp, 0.5s against 4.0s; and on 100,000 amplicons paired-end at
250bp on MSv3, 1.5s against 9.2s.

Apart from the indel fix described below, none of this changes the reads: given
the same seed and `-ir 0 -dr 0 -ir2 0 -dr2 0`, the patched and the stock
art_illumina produce byte-identical FASTQ.

The patch also gives art_illumina gzip and pipe support, which gargammel.pl now
uses to write the reads compressed in one pass instead of shelling out to gzip
afterwards. The same applies at the other end of that stage: because the
patched art reads a gzipped reference, adptSim writes the amplicons straight to
`<prefix>_a.fa.gz` (any `-arts`/`-artp` destination ending in `.gz` is
compressed as it is written) and art reads that file as it is. The amplicons
never exist uncompressed and are no longer gzipped in a pass of their own after
art has finished.

gargammel.pl asks for compression level 4 rather than gzip's
default of 6: on a simulated FASTQ, level 6 costs about five times the
compression time of level 4 for 7% off the file size (5.98s and 11.1MB against
1.22s and 12.0MB, on 200,000 reads here). The reads themselves are unaffected,
only the size of the container; change the `-gzl 4` in gargammel.pl if you want
the smaller files back.

* gzipped FASTA input is detected and decompressed with no flag;
* `-gz` (and `-gzl [1-9]` for the level) compresses the FASTQ/ALN/SAM output;
* `-i` and the new `--fq1`, `--fq2`, `--aln1`, `--aln2` and `--samFile` options
  accept `-`, `/dev/stdin`, `/dev/stdout`, `/dev/fd/N` and `fd:N`, so ART can be
  put in the middle of a pipeline:

`gzip -cd amplicons.fa.gz | art_illumina -ss HS25 -amp -na -p -l 75 -c 1 -i - --fq1 r1.fq.gz --fq2 r2.fq.gz`

When the reads go to stdout the run summary is written to stderr instead, so
the two never mix. Building the ALN or SAM header needs a second pass over the
reference, which is impossible on a pipe; ART says so rather than writing a
headerless file, so use `-na` without `-sam` when reading from stdin.



Sequencing indels:
-------------------------------------------------------------------------------------

Besides substituting bases, art_illumina also inserts and deletes them, at a low
per-base rate, so that roughly one read in eighty at 75bp carries an indel. What
those rates are and where they come from is set out below, along with when you
might want to override them; `gargammel.pl --noindel` switches indels off
entirely.

Note that gargammel produced almost no indels before version 1.1.5, because of a
bug in ART:

art_illumina's `-ir`/`-dr` insertion and deletion rates were not being realized
in ART 2.5.8: the table in `set_rate()` (`seqRead.h`) was built from P(X > i)
starting at i=1, so index i held P(X >= i+2) and placing one indel was gated on
the probability of needing two. At the default rates, 0.003% of 75bp reads
carried an indel instead of the ~1.5% the parameters imply. The patch starts
that table at i=0, and the measured rate is now 1.277% on gargammel amplicons,
matching an independent reimplementation of ART (`art_modern` 1.5.1, 1.275% on
the same input).

This means reads simulated with gargammel 1.1.5 differ from those of 1.1.4 even
at the same seed: about one read in eighty now carries a sequencing indel.
`gargammel.pl --noindel` turns them off, which reproduces the
substitution-only reads earlier versions produced:

`gargammel.pl -n 1000000 --comp 0,0,1 -l 40 -rl 75 --noindel -o data/simulation data/`

### Where the indel rates come from, and when to change them

The rates are art_illumina's own defaults. They are per base, not per read, and
the reverse read is given roughly twice the rate of the forward one:

| | insertion | deletion |
|---|---|---|
| forward read (`-ir`, `-dr`) | 9e-5 | 1.1e-4 |
| reverse read (`-ir2`, `-dr2`) | 1.5e-4 | 2.3e-4 |

Since they apply per base, the chance that a read carries at least one indel
grows with the read length you ask for with `-rl`:

| read length | forward | reverse |
|---|---|---|
| 35bp | 0.70% | 1.32% |
| 75bp | 1.49% | 2.81% |
| 100bp | 1.98% | 3.73% |
| 150bp | 2.96% | 5.54% |
| 250bp | 4.88% | 9.06% |

Simulated runs land close to that: counting gapped alignments in art_illumina's
own ALN output gives 1.28% of forward reads at 75bp and 3.06% at 150bp, against
the 1.49% and 2.96% above. The small discrepancies are sampling, plus the fact
that art draws once per entry of its indel table rather than inverting the
distribution in one go.

The only published account of where these four numbers come from is one
sentence in the ART paper: "the built-in insertion and deletion error rates
were derived from 35 bp reads aligned with our modified ACANA alignment tool".
The values themselves do not appear in the paper, and they are hardcoded in
art_illumina with no comment or citation. Three things follow, and they matter
if the indel rate is something your analysis is sensitive to:

* they were estimated on **35bp reads on the Illumina chemistry of around
  2010**, so they are not a description of any current platform;
* they are the same **whichever platform you select with `-ss`**. This is
  unlike the substitution errors, which are drawn from the empirical quality
  profile of the platform you chose, and so do track it;
* the per-read numbers above follow from the per-base rate and nothing else.
  A 250bp run does not have a measured indel rate of 4.9%; that is simply what
  9e-5 and 1.1e-4 per base come to over 250 bases.

Treat them as a plausible default rather than a calibrated one. If you have an
estimate for the data you are emulating — mapping your real reads and counting
gapped alignments will give you one — pass it to art_illumina directly with
`-ir`, `-dr`, `-ir2` and `-dr2`, or switch indels off with
`gargammel.pl --noindel` if they would confound what you are measuring.

    Huang, Weichun, et al. "ART: a next-generation sequencing read simulator."
    Bioinformatics 28.4 (2012): 593-594.


Specifying damage/deamination:
-------------------------------------------------------------------------------------

If you use gargammel.pl or deamSim, you can speficiy deamination/damage using either:

1. Use Briggs model parametes (see Briggs, Adrian W., et al. "Patterns of damage in genomic DNA sequences from a Neandertal." Proceedings of the National Academy of Sciences 104.37 (2007): 14616-14621.)

    -damage v,l,d,s

The four parameters and what they do to a molecule are described in the section
on the Briggs model below.

2. Use a misincorporation matrix computed by mapDamage (https://ginolhac.github.io/mapDamage). This matrix is in the results directory created by mapDamage and is called "misincorporation.txt". There are 2 examples of such files:

    examplesMapDamage/results_LaBrana/misincorporation.txt
    examplesMapDamage/results_Ust_Ishim/misincorporation.txt

The first is from a  double-stranded library and the second a single-stranded one. To use either, you can use the wrapper script or deamSim as such:

    -mapdamage examplesMapDamage/results_LaBrana/misincorporation.txt double
    -mapdamage examplesMapDamage/results_Ust_Ishim/misincorporation.txt single

We suggest that you run mapDamage on the empirical data that you are trying to emulate and use the resulting misincorporation.txt file.

3. Specify a matrix of deamination rates, we use the following format, the first line is the header:

    	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
    	pos	rate_{A->C}	rate_{A->G}	rate_{A->T}	rate_{C->A}	rate_{C->G}	rate_{C->T}	rate_{G->A}	rate_{G->C}	rate_{G->T}	rate_{T->A}	rate_{T->C}	rate_{T->G}

The pos. is the position 0,1... after the fragment beginning/end. The rate is specified using the following: estimate  [estimate_low estimate_high]. For example, 0.3 [0.2 0.4] means that the rate of deamination is 0.3 or 30%.

example of a format:

	A->C	A->G	A->T	C->A	C->G	C->T	G->A	G->C	G->T	T->A	T->C	T->G
	0	1.853e-3 [1.726e-3..1.989e-3]	4.064e-3 [3.875e-3..4.263e-3]	3.269e-3 [3.099e-3..3.448e-3]	6.661e-3 [6.254e-3..7.094e-3] 3.057e-3 [2.785e-3..3.355e-3] 8.004e-2 [7.865e-2..8.145e-2] 1.236e-2 [    1.183e-2..1.292e-2] 4.131e-3 [3.828e-3..4.459e-3] 6.703e-3 [6.314e-3..7.116e-3] 3.845e-3 [3.624e-3..4.079e-3] 4.581e-3 [4.339e-3..4.836e-3] 2.169e-3 [2.005e-3..2.347e-3]
	1	1.986e-3 [1.849e-3..2.134e-3]	4.273e-3 [4.070e-3..4.487e-3]	3.030e-3 [2.859e-3..3.211e-3]	5.357e-3 [5.001e-3..5.738e-3] 3.188e-3 [2.916e-3..3.485e-3] 1.427e-2 [1.369e-2..1.488e-2] 9.514e-3 [    9.075e-3..9.974e-3]	3.316e-3 [3.061e-3..3.593e-3] 5.061e-3 [4.743e-3..5.400e-3] 3.421e-3 [3.216e-3..3.639e-3] 4.865e-3 [4.620e-3..5.124e-3]	2.201e-3 [2.038e-3..2.377e-3]

This follows the output of https://bitbucket.org/ustenzel/damage-patterns.git

4. You can use one of the precalculated rates of deamination in src/matrices/. There is a damage from single-strand and a double-strand libraries from the following study: 

    Lazaridis, Iosif, et al. "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513.7518 (2014): 409-413.

See the methylation question for adding different rates of deamination for methylated/unmethylated cytosine.


The Briggs model:
-------------------------------------------------------------------------------------

    -damage v,l,d,s

This is the model of:

    Briggs, Adrian W., et al. "Patterns of damage in genomic DNA sequences from a Neandertal." Proceedings of the National Academy of Sciences 104.37 (2007): 14616-14621.

Rather than giving a rate of C->T per position, it describes the physical state
of an ancient molecule and lets the damage pattern follow from it. Four
parameters are needed:

| parameter | meaning |
|---|---|
| `v` | probability that the molecule carries a nick |
| `l` | geometric parameter for the length of the single-stranded overhangs |
| `d` | probability that a cytosine in a double-stranded region is deaminated |
| `s` | probability that a cytosine in a single-stranded region is deaminated |

The maximum likelihood estimates of Briggs et al. for their Neandertal data are
`-damage 0.024,0.36,0.0097,0.68`, which is a reasonable starting point for a
double-stranded library.

### What the model says about one molecule

An ancient molecule is double-stranded in the middle and frayed at the ends,
where one strand extends past the other. Cytosine deaminates to uracil far more
readily when it is not paired, which is why the damage concentrates at the
ends: the single-stranded overhangs deaminate at `s`, the double-stranded
interior only at the much lower `d`. The uracils are read as thymine, so a
deaminated C appears as T.

The length of each overhang is drawn from a geometric distribution with
parameter `l`: `l` is the probability of stopping at each base, so a larger `l`
means shorter overhangs, and an overhang of length zero has probability `l`.
An overhang is as likely to leave the 5' strand protruding as the 3' one, and
the blunt-end repair step of the library preparation fills in the 3' overhangs
while leaving the 5' ones. Each end of the molecule therefore shows a
single-stranded region only half of the time.

Which substitution you see depends on which end you are at. A 5' overhang is
sequenced as it is, so its deaminated cytosines show up as C->T. A 3' overhang
is read on the complementary strand, where the same deamination appears as
G->A. That asymmetry — C->T at the 5' end, G->A at the 3' end — is the
signature of a double-stranded library, and it is what the model produces
without being told to.

Nicks are the second ingredient. A molecule with a nick, which happens with
probability `v`, is copied from the nick onwards along the other strand, so
everything downstream of the nick is read in the opposite orientation: the
double-stranded interior after the nick shows G->A instead of C->T. The nick
sits uniformly within the double-stranded region. A nick inside an overhang
would merely shorten it and is not observable, and Briggs et al. note that
molecules with a second nick on the opposite strand are lost during the repair,
which "causes the distribution of first nicks in the sequenced fragments to be
uniform rather than geometric".

Molecules so short that the two overhangs meet are single-stranded over their
whole length and deaminate at `s` everywhere.

### The rate you should expect at the first base

For the estimates above, the 5'-most base is inside a 5' overhang with
probability 0.5*(1-l) = 0.32 and deaminates there with probability s = 0.68;
the rest of the time it is double-stranded and deaminates with probability
d = 0.0097, which works out at about 22%. Simulated fragments come out at
21-22% C->T at the first position and the same G->A at the last one, which is
the ~21% reported by Briggs et al. The rate then falls off geometrically with
the distance from the end, by a factor of about (1-l) per base:

    position from 5' end     1      2      3      4      5      6
    C->T                     0.215  0.143  0.099  0.067  0.043  0.033

For comparison, the pre-1.1.5 behavior kept by `-damagelegacy` starts at 0.440
for the same parameters.

### Note on versions before 1.1.5

Two aspects of this model were implemented incorrectly up to gargammel 1.1.4 and were corrected afterwards:

  * The factor of one half above was missing: gargammel drew a geometric overhang at both ends unconditionally, as though the blunt-end repair never removed one. This doubled the rate of C->T and G->A at the terminal positions. With the parameters of Briggs et al. the 5'-most C->T rate came out at 44% instead of the 21% reported in that paper.

  * The nick was placed geometrically rather than uniformly, and it was looked for across the single-stranded overhangs as well, where a nick is not observable. This shifted the C->T/G->A crossover in the interior of the molecule towards the 5' end. The effect is bounded by `d` and is therefore much smaller than the one above.

If you need to reproduce a simulation made with an earlier version, use `-damagelegacy` instead of `-damage` (or pass `--damagelegacy` to gargammel.pl, which switches every `-damage*` option at once). It takes the same four parameters and restores the previous behavior exactly.


Can I specify different rates of misincorporation due to deamination for the endogenous/bacterial/human contaminant sources?
-------------------------------------------------------------------------------------

Yes, please refer to the options of the wrapper script gargammel.pl


Is it possible to specify different rates of deamination for methylated and unmethylated bases?
-------------------------------------------------------------------------------------

Yes. In the endogenous genome, specify methylated cytosine as 'c' (lowercase c) and unmethylated cytosine as 'C' (uppercase C). You can specify multiple cells using the following directory structure:

    input/
    input/endo
    input/endo/C0
    input/endo/C0/chr20_0_split1.fa
    input/endo/C0/chr20_0_split1.fa.fai
    input/endo/C0/chr20_0_split2.fa
    input/endo/C0/chr20_0_split2.fa.fai
    input/endo/C1
    input/endo/C1/chr20_1_split1.fa
    input/endo/C1/chr20_1_split1.fa.fai
    input/endo/C1/chr20_1_split2.fa
    input/endo/C1/chr20_1_split2.fa.fai
    input/endo/C2
    input/endo/C2/chr20_2_split1.fa
    input/endo/C2/chr20_2_split1.fa.fai
    input/endo/C2/chr20_2_split2.fa
    input/endo/C2/chr20_2_split2.fa.fai

Where C0 reprensents the first cell, C1 the second and so forth. A lower case C 'c' is a methylated C and an uppercase 'C' is a an unmethylated C. To create these files from a reference and a methylation map, please see the script src/addMethyl.pl which needs to be modified (hardcoded paths).

Methylated and unmethylated cytosines on the - strand can be specified using 'g' and 'G'. Once this is done, you can specify the option:  --methyl for gargammel.pl.  When using --methyl, you can then specify different matrix files for rates of deamination for nonmethylated and methylated cytosines:

    -matfilenonmeth    [matrix file prefix] Read the matrix file of substitutions for non-methylated Cs
    -matfilemeth       [matrix file prefix] Read the matrix file of substitutions for methylated Cs


How can I get an ancient DNA composition profile for gargammel?
-------------------------------------------------------------------------------------

By composition we mean the base frequency at the breaks. You could generate it manually, the format is as follows:

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

	Jonsson, Hakon, et al. "mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters." Bioinformatics (2013): btt193.

It is normally called "dnacomp.txt" in the output directory, you can filter a single chromosome (in this case 21) using this command:

	grep "^21\|^#\|^Chr"  /path to mapDamage output/results_[sample name]/dnacomp.txt >  dnacomp.txt


How can I specify the size distribution?
-------------------------------------------------------------------------------------

Ancient DNA molecules tend to be fragmented and can be very short but tend to have a specific shape. Both for the wrapper script (gargammel.pl) and the fragment simulation program (fragSim), there are are 4 ways to specify the :

1) Specify a fixed length using -l 
------------

2) Open a file containing the size distribution using -s, one empirical fragment length per line eg:
------------
~~~~
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
~~~~
------------
3) Open a file containing the size frequencies using -f in the format "size[TAB]freq" eg:
------------
~~~~
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
~~~~

------------
4) Specify the size distribution using parameters from a log-normal distribution, using options --loc and --scale.
------------


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

For the input/bact/ directory which represent the microbial contamination, gargammel needs a set of fasta files that represent the different microbes. **Each file corresponds to exactly one microbial species.** Each fasta file must contain the genome of the microbial species, multiple scaffolds and plasmids are allowed. Each fasta file must also be faidx indexed. This directory must also contain a file called "list".  This file contains the list of every fasta files in that directory along with their relative abundance in the desired bacterial contamination. For example:

    bacteria1.fa	0.5
    bacteria2.fa	0.3
    bacteria3.fa	0.2

The abundance will be printed on the console when the program is launched. Some users have reported discrepancies between the original bacterial abundance and the printed one. Make sure that they are equal and that the bacterial abundance file uses UNIX carriage returns (use dos2unix or mac2unix to transform from DOS/MAC to Unix format).

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


This will download the necessary files from NCBI to create a database suitable for gargammel to simulate microbial species in the exampleBacteriaDB/fasta and run samtools faidx on each file. You need standard UNIX utilities such as awk/sed/python/curl/wget/gzip to be installed as well as samtools. Please move the fasta/ directory produced (exampleBacteriaDB/fasta in the example above) to the input/bact/. The file named "exampleBacteriaDB/fastafasta/list" is the list of bacterial species along with their abundance. Another file, "exampleBacteriaDB/Microbial_ID.log" details the strain/ID and ftp link used. retrieveFromMetabit uses GNU parallel (see O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.), please make sure that it is installed.

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

FAQ and issues 
-------------------------------------------------------------------------------------

* I am getting:

    ./art_illumina: error while loading shared libraries: libgsl.so.0: cannot open shared object file: No such file or directory

Make sure you have libgsl installed and create a symbolic link:

     sudo ln -s  /usr/lib/x86_64-linux-gnu/libgsl.so.23.0.0   /usr/lib/libgsl.so.0 

* Can I know where on the genome of origin the fragment was sampled?

    Yes! fragSim, which generates the original fragments uses the following format: [CHROMOSOME NAME]:[STRAND]:[START]:[END]:[LENGTH]
    The overall wrapper (gargammel.pl) will add e1_ to endogenous fragments from the first reference and e2_ to endogenous fragments from the second reference. It will add c_X to the fragments from present-day human contaminants where X is the # of the genome and b_X to the fragments from bacterial contaminants where X is the # of the genome.

* Can gargammel simulate indels?

   Yes and no. gargammel does not currently insert indels as a result of sequencing errors. However, if you add indels in your input genome, it will handle them without any problems. 
