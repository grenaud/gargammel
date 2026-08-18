gargammel release notes
=====================================================================================

1.1.5
-------------------------------------------------------------------------------------

This release corrects the implementation of the Briggs model and a long-standing
indel bug in the bundled read simulator, makes a whole simulation reproducible
from a single seed, and ships a patched art_illumina that is several times
faster and can read and write compressed streams.

**This release changes the reads you get.** Two corrections are responsible,
and both are deliberate:

* the damage at the terminal positions was too high by a factor of two under
  `-damage`, and is now roughly half what it was for the same parameters.
  Earlier simulations can still be reproduced exactly with `-damagelegacy`;
* sequencing indels were essentially absent — 0.003% of reads instead of the
  ~1.5% the ART parameters imply — and now appear at the rate asked for.

Both are switchable: `-damagelegacy` restores the old damage model and
`--noindel` turns the sequencing indels back off, so passing the two together
puts 1.1.5 back on the code path 1.1.4 used, and its output is comparable with
what you measured before. Note that this is not bit-for-bit reproduction of a
particular old run: nothing before 1.1.5 was seedable, so the random stream
behind an existing pre-1.1.5 dataset cannot be recovered.


### The Briggs model is now implemented as described in the paper

`-damage v,l,d,s` simulates the model of Briggs et al. (2007). Two aspects of
it were implemented incorrectly up to 1.1.4:

* **Overhangs were drawn at both ends unconditionally.** Briggs et al. assume
  an overhang is equally likely to leave the 5' or the 3' strand protruding,
  and that the blunt-end repair step of the library preparation fills in the 3'
  overhangs while leaving the 5' ones, so each end of a molecule carries an
  observable single-stranded region only half of the time. gargammel drew a
  geometric overhang at both ends every time, which doubled the rate of C->T
  and G->A at the terminal positions. With the maximum likelihood estimates of
  Briggs et al. (`-damage 0.024,0.36,0.0097,0.68`) the 5'-most C->T rate came
  out at 44% instead of the ~21% reported in that paper.

* **The nick was placed geometrically, and looked for inside the overhangs.**
  Briggs et al. note that molecules carrying a second nick on the opposite
  strand are lost during the nick repair, which "causes the distribution of
  first nicks in the sequenced fragments to be uniform rather than geometric".
  A nick falling inside a single-stranded overhang merely shortens it and is
  not observable. Placing it geometrically and scanning across the overhangs
  shifted the C->T/G->A crossover in the interior of the molecule towards the
  5' end. The size of this effect is bounded by `d`, so it is much smaller than
  the one above.

Simulated with the estimates of Briggs et al., the C->T rate by distance from
the 5' end is now:

    position from 5' end     1      2      3      4      5      6
    1.1.5                    0.215  0.143  0.099  0.067  0.043  0.033
    1.1.4 and earlier        0.440  0.277  0.186  0.119  0.078  0.056

The first column is the ~21% of the paper. The 3' G->A rates are the mirror
image of these, as they should be.

`doc/briggs_model_note.tex` derives both corrections step by step, and the
README now has a section describing the model itself: what the four parameters
mean, why the damage concentrates at the ends, and why a double-stranded
library shows C->T at the 5' end and G->A at the 3' end.

**Reproducing older simulations.** `deamSim -damagelegacy v,l,d,s` takes the
same four parameters and restores the previous behavior exactly.
`gargammel.pl --damagelegacy` switches every `-damage*` option at once. This
exists only for reproducing published runs; use `-damage` for new work.


### Reproducible simulations from a single seed

`gargammel.pl --seed [int]` makes a whole run byte-reproducible: two runs with
the same seed and the same arguments produce identical output. Without it,
every run is seeded from the clock, exactly as before.

The wrapper seeds its own draws — it splits the requested number of fragments
over the input genomes at random — and hands each of fragSim, deamSim, adptSim
and art a *separate* seed derived from the master one, so the sub-programs do
not share a random stream. Each of them also takes its own `--seed` (`-rs` for
art) if you drive it directly.

Making this hold required two changes:

* **adptSim now takes `--seed`.** It pads fragments shorter than the read
  length with random bases, and that padding came from a clock-seeded
  generator, which was on its own enough to make a run irreproducible.

* **art_illumina now carries its own random number generator** rather than
  calling the one in the C library. It reproduces glibc's `rand()` exactly, so
  runs match results from earlier versions on Linux, and art_illumina on its
  own now gives the same reads from the same seed on Linux and on macOS.

That cross-platform guarantee covers art_illumina only. fragSim, deamSim and
adptSim still draw through `rand`, `drand48` and `std::default_random_engine`,
all of which differ between implementations, so a whole gargammel run with a
fixed seed reproduces on the same machine and toolchain, not across platforms.


### A faster art_illumina, with gzip and pipe support

ART is the slowest stage of the pipeline. The copy gargammel builds is now
patched — see `patches/art_illumina_gargammel.patch`, applied by the Makefile
when ART is unpacked — and runs several times faster than the stock version.
Measured against a stock art_illumina built from the same tarball:

| run | stock | patched |
|---|---|---|
| 200k amplicons, paired-end, 75bp, HS25 | 7.3s | 0.9s |
| 200k amplicons, single-end, 75bp, HS25 | 4.0s | 0.5s |
| 100k amplicons, paired-end, 250bp, MSv3 | 9.2s | 1.5s |

Every part of the patch except the indel fix below leaves the reads untouched:
run the patched and the stock version with the same seed and `-ir 0 -dr 0 -ir2 0
-dr2 0` and the FASTQ files are byte-identical.

The patch also gives art_illumina compressed and non-seekable I/O:

* gzipped FASTA input is detected and decompressed with no flag;
* `-gz` (and `-gzl [1-9]` for the level) compresses the FASTQ/ALN/SAM output;
* `-i` and the new `--fq1`, `--fq2`, `--aln1`, `--aln2` and `--samFile` options
  accept `-`, `/dev/stdin`, `/dev/stdout`, `/dev/fd/N` and `fd:N`, so ART can
  be put in the middle of a pipeline. When the reads go to stdout the run
  summary goes to stderr instead, so the two never mix.

Building an ALN or SAM header needs a second pass over the reference, which is
impossible on a pipe; ART now says so rather than writing a headerless file, so
use `-na` without `-sam` when reading from stdin.


### Sequencing indels are now simulated at the rate ART is asked for

art_illumina takes an insertion and a deletion rate per base (`-ir`, `-dr`, and
`-ir2`/`-dr2` for the second read; the defaults are 9e-5 and 1.1e-4). Those
rates were not being realized: at 75bp, 0.003% of reads carried an indel rather
than the ~1.5% the parameters imply, roughly 250x too few. In practice
gargammel simulated substitution errors and nothing else.

The cause is in `set_rate()` in ART's `seqRead.h`, which builds the table
`get_indel()` uses to decide how many indels a read gets. Reading that table
from the top down, index `i` places `i+1` indels, so index `i` has to hold
P(X >= i+1). The loop filling it started at `i=1` and stored
`gsl_cdf_binomial_Q(i, p, read_len)` = P(X > i), so index `i` actually held
P(X >= i+2): placing a *single* indel was gated on the probability of needing
*two*. The `i=0` term was computed, used to seed the loop's cutoff accumulator,
then overwritten without ever being stored. The patch starts the table at
`i=0`.

Measured on 200,000 amplicons at 75bp, reads carrying an indel:

| | rate |
|---|---|
| 1.1.4 and earlier | 0.003% |
| 1.1.5 | 1.277% |
| nominal, from `-ir`/`-dr` | ~1.5% |

As an independent check, `art_modern` 1.5.1 — a separate reimplementation of
ART by YU Zhejian, which does not have this bug — gives 1.275% on the same
input.

The rates are art_illumina's own, per base: 9e-5 insertions and 1.1e-4
deletions on the forward read, 1.5e-4 and 2.3e-4 on the reverse. The ART paper
says only that they "were derived from 35 bp reads aligned with our modified
ACANA alignment tool" and does not publish the values; they are hardcoded with
no comment. They are therefore not tied to the platform chosen with `-ss` — the
substitution errors are, since those come from the empirical quality profile —
and the chance of a read carrying one rises with `-rl`, the rate being per
base. Use `-ir`/`-dr`/`-ir2`/`-dr2` on art_illumina directly if you have an
estimate for your own data.

Two consequences worth knowing:

* **Reads from 1.1.5 are not byte-identical to reads from 1.1.4**, even with
  the same seed and the same damage settings, because roughly one read in
  eighty now carries an indel that it did not carry before. Everything else
  about the read is unchanged. The new `gargammel.pl --noindel` zeroes all four
  rates and reproduces the old substitution-only reads.
* Benchmarks built on gargammel output may get slightly harder. Aligners had
  been seeing reads that differed from the reference by substitutions only,
  which is easier than real data; mapping rates measured against earlier
  versions will not be directly comparable.

This is a defect in ART 2.5.8 itself, not something gargammel introduced, so it
affects other tools that wrap this version of art_illumina.


### A --noindel flag

`gargammel.pl --noindel` passes `-ir 0 -dr 0 -ir2 0 -dr2 0` to art_illumina, so
the only sequencing error is substitution. This is the escape hatch for the
indel change above: it gives the substitution-only reads earlier versions
produced, and is also useful when indels would confound whatever you are
measuring.


### No more separate gzip passes

gargammel.pl used to write the reads and the amplicons uncompressed and then
shell out to `gzip` once ART had finished. Both files are now written
compressed in one pass:

* the reads are compressed by art itself, at level 4 rather than gzip's default
  of 6 — on a simulated FASTQ, level 6 costs about five times the compression
  time of level 4 for 7% off the file size. Change the `-gzl 4` in
  `gargammel.pl` if you prefer the smaller files;
* the amplicons are written straight to `<prefix>_a.fa.gz` by adptSim and read
  as such by the patched art, so they never exist uncompressed. Any
  `-arts`/`-artp` destination ending in `.gz` is now compressed as it is
  written; any other name is written plain, as before.

The output file names are unchanged.


### Static binaries

    make static

builds everything the ordinary way, then relinks the six programs in `src/` and
`art_illumina` against the static libc, libstdc++, libz and libgsl, so the
binaries can be copied to a machine that does not have them. `gargammel.pl`
itself is a Perl script and still needs perl. A static build produces the same
reads as an ordinary one, seed for seed.

This needs the static system libraries, which on most distributions are a
separate package from the headers: on Debian/Ubuntu they come with `libc6-dev`,
`zlib1g-dev` and `libgsl-dev`. macOS does not ship a static libc, so use the
ordinary build there. Run `make clean` before going back to dynamic binaries —
a plain `make` will leave the static ones in place, since they are newer than
the sources.


### Documentation

* A new section shows how to run the whole simulation as **one command**,
  chaining fragSim, deamSim, adptSim and art_illumina through pipes with no
  intermediate files, for the common case of drawing from a single genome. It
  covers both paired-end and single-end, and shows how to drop a `tee` between
  any two stages to keep the intermediate files — the undamaged fragments, the
  damaged fragments, the amplicons — without interrupting the flow.
* A new section describes **the Briggs model** itself rather than only its
  options.
* New sections on reproducible simulations, on the speed and I/O of the read
  simulator, and on static binaries.


### Tests

`make test` now runs 299 checks over ten groups, up from 247 over nine. The new
`art` group covers the patched art_illumina — seeding, gzip in and out, the
pipe forms, and that a SAM run off a pipe is refused rather than silently
producing a headerless file. The `deamSim` group checks the Briggs model
against the rates published by Briggs et al. and that `-damagelegacy` still
reproduces the old ones; the `gargammel` group checks that a fixed seed
reproduces a run exactly, that every sub-program gets a distinct seed, and that
neither the reads nor the amplicons are gzipped in a pass of their own.

As before, the suite is deterministic, builds its own inputs, writes only to a
temporary directory and needs no data beyond this repository.


### Upgrading

* Simulations that use `-damage` will differ from 1.1.4: the terminal damage
  rates are now about half of what they were. This is the point of the release.
  Use `-damagelegacy` with the same parameters to reproduce an older run.
* Simulations that specify damage with `-mapdamage`, `-matfile` or `-profile`
  are unaffected *by the damage change*, but every simulation is affected by
  the indel fix, since that happens downstream in art_illumina.
* Roughly one read in eighty now carries a sequencing indel that it would not
  have carried before. If you are comparing mapping rates or error profiles
  against numbers measured with an earlier version, regenerate them.
* `<prefix>_a.fa.gz` is now written by adptSim rather than by a `gzip` pass, and
  is present while art is running rather than appearing at the end. Its name
  and contents are unchanged.
* Anything that called the bundled `art_illumina` directly still works. The
  patch only adds options, and the only change to its output is the indel fix,
  which `-ir 0 -dr 0 -ir2 0 -dr2 0` turns off.


### References

Briggs, Adrian W., et al. "Patterns of damage in genomic DNA sequences from a
Neandertal." *Proceedings of the National Academy of Sciences* 104.37 (2007):
14616-14621.
