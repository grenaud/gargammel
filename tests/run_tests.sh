#!/usr/bin/env bash
#
# Test suite for the gargammel subcomponents.
#
# Each of the compiled tools (fragSim, deamSim, adptSim, fasta2fastas,
# mapDamage2prof, damage_patterns2prof) and the gargammel.pl wrapper is run on
# small synthetic inputs generated on the fly and its output is checked.
#
# Usage:  tests/run_tests.sh [options]
#
#   --only [pattern]   Only run the test groups matching [pattern]. Several
#                      names may be given, separated by | or by spaces, and
#                      each may be a shell wildcard, e.g.
#                      --only 'fragSim|deamSim'   --only '*Sim'
#   --list             List the test groups and exit
#   --keep             Do not delete the temporary working directory
#   -v, --verbose      Echo every command that is run
#   -h, --help         This message
#
# Returns 0 if every test passed, 1 otherwise.
#

TESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GARGDIR="$( dirname "$TESTDIR" )"
SRCDIR="$GARGDIR/src"

TESTGROUPS="usage fasta2fastas fragSim deamSim adptSim mapDamage2prof damage_patterns2prof pipeline gargammel"

ONLY=""
KEEP=""
VERBOSE=""

while [ $# -gt 0 ]; do
    case "$1" in
	--only)    ONLY="$2"; shift 2 ;;
	--list)    echo $TESTGROUPS | tr ' ' '\n'; exit 0 ;;
	--keep)    KEEP=1; shift ;;
	-v|--verbose) VERBOSE=1; shift ;;
	-h|--help) sed -n '2,20p' "$0" | sed 's/^# \?//'; exit 0 ;;
	*)         echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [ -t 1 ]; then
    CGREEN=$'\033[32m'; CRED=$'\033[31m'; CYELLOW=$'\033[33m'; CBOLD=$'\033[1m'; COFF=$'\033[0m'
else
    CGREEN=""; CRED=""; CYELLOW=""; CBOLD=""; COFF=""
fi

NPASS=0
NFAIL=0
NSKIP=0
FAILED=""

W=$(mktemp -d "${TMPDIR:-/tmp}/gargammel-test.XXXXXX") || exit 1
if [ -z "$KEEP" ]; then
    trap 'rm -rf "$W"' EXIT
fi

pass() { NPASS=$((NPASS+1)); printf '  %sok%s   %s\n' "$CGREEN" "$COFF" "$1"; }
fail() {
    NFAIL=$((NFAIL+1))
    FAILED="$FAILED
    $1"
    printf '  %sFAIL%s %s\n' "$CRED" "$COFF" "$1"
    [ -n "$2" ] && printf '         %s\n' "$2"
    return 0
}
skip() { NSKIP=$((NSKIP+1)); printf '  %sskip%s %s (%s)\n' "$CYELLOW" "$COFF" "$1" "$2"; }
group() { printf '\n%s== %s%s\n' "$CBOLD" "$1" "$COFF"; }

# eq [description] [expected] [actual]
eq() {
    if [ "x$2" = "x$3" ]; then pass "$1"; else fail "$1" "expected [$2] but got [$3]"; fi
}

# ne [description] [not expected] [actual]
ne() {
    if [ "x$2" != "x$3" ]; then pass "$1"; else fail "$1" "did not expect [$3]"; fi
}

# ok [description] [return code] [extra info on failure]
ok() {
    if [ "$2" -eq 0 ]; then pass "$1"; else fail "$1" "${3:-returned $2}"; fi
}

# run [logname] [cmd...]  -- runs a command, stdout in $W/[logname].out,
# stderr in $W/[logname].err, sets $RC.  Nothing here should take minutes, so
# a command that does is killed rather than left to stall the whole suite.
TIMEOUT=300
run() {
    local name="$1"; shift
    [ -n "$VERBOSE" ] && printf '     $ %s\n' "$*" >&2
    if command -v timeout >/dev/null 2>&1; then
	timeout "$TIMEOUT" "$@" >"$W/$name.out" 2>"$W/$name.err"
    else
	"$@" >"$W/$name.out" 2>"$W/$name.err"
    fi
    RC=$?
    if [ "$RC" -eq 124 ]; then
	echo "timed out after ${TIMEOUT}s" >> "$W/$name.err"
    fi
    return $RC
}

# tail of the stderr of the last run, for failure messages
lasterr() { tail -n 2 "$W/$1.err" 2>/dev/null | tr '\n' ' '; }

# does the group name match what --only asked for? The pattern may hold
# several names separated by | or by spaces, each of them a shell wildcard
wants() {
    [ -z "$ONLY" ] && return 0
    local pat
    for pat in ${ONLY//|/ }; do
	case "$1" in $pat) return 0 ;; esac
    done
    return 1
}

#####################################################################
# fasta helpers.  Every tool here writes one sequence per line, so
# a record is a defline followed by a single sequence line.
#####################################################################

# cat a fasta whether it is gzipped or not
fa_cat() {
    case "$1" in
	*.gz) gzip -cd "$1" ;;
	*)    cat "$1" ;;
    esac
}

fa_count() { fa_cat "$1" | grep -c '^>' ; }
fa_seqs()  { fa_cat "$1" | grep -v '^>' ; }
fa_lens()  { fa_seqs "$1" | awk '{print length($0)}' ; }
# the sorted set of distinct sequence lengths, comma separated
fa_lenset() { fa_lens "$1" | sort -n -u | tr '\n' ',' ; }

# is the file a valid gzip stream?
is_gzip() { gzip -t "$1" >/dev/null 2>&1 ; }

# is the file a BAM file? (BGZF is valid gzip and starts with the BAM magic)
is_bam() {
    is_gzip "$1" || return 1
    [ "$(gzip -cd "$1" 2>/dev/null | head -c 4)" = "BAM$(printf '\1')" ]
}

#####################################################################
# Fixtures
#####################################################################

# A deterministic pseudo-random DNA sequence, so that every run of the
# suite works on exactly the same reference. MINSTD generator, the
# arithmetic stays well inside the 53 bits of an awk double.
make_fasta() {
    local out="$1" name="$2" len="$3" seed="$4"
    awk -v name="$name" -v len="$len" -v x="$seed" 'BEGIN{
	printf ">%s\n", name;
	line="";
	for(i=0;i<len;i++){
	    x=(16807*x)%2147483647;
	    line=line substr("ACGT", int(x/1024)%4+1, 1);
	    if(length(line)==60){ print line; line=""; }
	}
	if(length(line)>0) print line;
    }' > "$out"
    # the .fai that fragSim mmaps; deterministic since we control the layout
    awk -v name="$name" 'BEGIN{ off=length(name)+2 } !/^>/{ n+=length($0) } END{
	print name "\t" n "\t" off "\t60\t61"
    }' "$out" > "$out.fai"
}

setup_fixtures() {
    REF="$W/ref.fa"
    make_fasta "$REF" testchr 20000 42

    REF2="$W/ref2.fa"
    make_fasta "$REF2" testchr2 20000 4242

    # one fragment length per line, for fragSim -s
    printf '30\n40\n50\n30\n40\n50\n' > "$W/sizes.size"

    # length[TAB]frequency, for fragSim -f
    printf '25\t0.5\n45\t0.5\n' > "$W/sizes.freq"

    # a fasta carrying every IUPAC ambiguity code fasta2fastas knows about
    printf '>iupacchr\nACGTRYSWKMACGTRYSWKM\n' > "$W/iupac.fa"

    MAPDAMAGE="$GARGDIR/examplesMapDamage/results_LaBrana/misincorporation.txt"
    MAPDAMAGE2="$GARGDIR/examplesMapDamage/results_Ust_Ishim/misincorporation.txt"
    DNACOMP="$SRCDIR/dnacomp.txt"

    # substitution matrices in the format of damage-patterns, i.e. a header
    # of 13 tab separated fields and then "value [low..high]" per cell.
    # In a damage-patterns 3' matrix the *last* row is the base next to the
    # end of the molecule, so the damage grows down the file there.
    make_dp_dat "$W/dp-5p.dat" 6 fwd
    make_dp_dat "$W/dp-3p.dat" 7 rev
    DPROWS=20
}

# a plausible looking damage-patterns matrix: column 6 is C>T and column 7 is
# G>A, the chosen one decays away from the end of the molecule, the rest is flat
make_dp_dat() {
    awk -v dcol="$2" -v dir="$3" 'BEGIN{
	n=20;
	printf "\tA>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";
	for(pos=0;pos<n;pos++){
	    printf "%d", pos;
	    d = (dir=="rev") ? (n-1-pos) : pos;
	    for(c=1;c<=12;c++){
		if(c==dcol) v=0.3*exp(-d/3.0); else v=0.002;
		printf "\t%.4e [%.4e..%.4e]", v, v*0.9, v*1.1;
	    }
	    printf "\n";
	}
    }' > "$1"
}

#####################################################################
# usage: every binary must be there and must document itself
#####################################################################

test_usage() {
    group "usage / smoke test"
    local prog
    for prog in fragSim deamSim adptSim fasta2fastas mapDamage2prof damage_patterns2prof; do
	if [ ! -x "$SRCDIR/$prog" ]; then
	    fail "$prog is built" "$SRCDIR/$prog is missing, run 'make' first"
	    continue
	fi
	pass "$prog is built"

	run "usage-$prog" "$SRCDIR/$prog"
	eq "$prog with no argument exits 1" 1 "$RC"

	run "help-$prog" "$SRCDIR/$prog" -h
	if grep -q -- "-h" "$W/help-$prog.out" || [ -s "$W/help-$prog.out" ]; then
	    pass "$prog -h prints a usage message"
	else
	    fail "$prog -h prints a usage message" "nothing on stdout"
	fi
    done
}

#####################################################################
# fasta2fastas: splits the IUPAC codes of a chromosome into two haplotypes
#####################################################################

test_fasta2fastas() {
    group "fasta2fastas"
    local prog="$SRCDIR/fasta2fastas"
    [ -x "$prog" ] || { skip "fasta2fastas" "not built"; return; }

    run f2f "$prog" "$W/iupac.fa" "$W/hap"
    ok "runs on a fasta with IUPAC codes" "$RC" "$(lasterr f2f)"

    local h1="$W/hap1.fa.gz" h2="$W/hap2.fa.gz"
    if [ ! -s "$h1" ] || [ ! -s "$h2" ]; then
	fail "writes hap1.fa.gz and hap2.fa.gz" "one of the two is missing or empty"
	return
    fi
    pass "writes hap1.fa.gz and hap2.fa.gz"

    if is_gzip "$h1" && is_gzip "$h2"; then pass "both haplotypes are gzipped"
    else fail "both haplotypes are gzipped"; fi

    eq "haplotype 1 keeps the defline" ">iupacchr" "$(fa_cat "$h1" | head -1)"
    eq "haplotype 2 keeps the defline" ">iupacchr" "$(fa_cat "$h2" | head -1)"

    local inseq s1 s2
    inseq=$(fa_seqs "$W/iupac.fa" | tr -d '\n')
    s1=$(fa_seqs "$h1" | tr -d '\n')
    s2=$(fa_seqs "$h2" | tr -d '\n')

    eq "haplotype 1 has the length of the input" "${#inseq}" "${#s1}"
    eq "haplotype 2 has the length of the input" "${#inseq}" "${#s2}"

    # at an unambiguous base both haplotypes must carry that base, at an
    # ambiguous one they must carry the two bases the code stands for, in
    # either order since fasta2fastas picks the phase at random
    local i c a b bad_fixed=0 bad_iupac=0 pair
    for (( i=0; i<${#inseq}; i++ )); do
	c=${inseq:$i:1}; a=${s1:$i:1}; b=${s2:$i:1}
	case "$c" in
	    R) pair="AG" ;; Y) pair="CT" ;; S) pair="CG" ;;
	    W) pair="AT" ;; K) pair="GT" ;; M) pair="AC" ;;
	    *) [ "$a" = "$c" ] && [ "$b" = "$c" ] || bad_fixed=$((bad_fixed+1)); continue ;;
	esac
	[ "$a$b" = "$pair" ] || [ "$b$a" = "$pair" ] || bad_iupac=$((bad_iupac+1))
    done
    eq "unambiguous bases are copied to both haplotypes" 0 "$bad_fixed"
    eq "IUPAC codes are resolved into their two bases" 0 "$bad_iupac"

    run f2f-args "$prog" "$W/iupac.fa"
    eq "rejects a wrong number of arguments" 1 "$RC"
}

#####################################################################
# fragSim: draws aDNA fragments from a chromosome
#####################################################################

test_fragSim() {
    group "fragSim"
    local prog="$SRCDIR/fragSim"
    [ -x "$prog" ] || { skip "fragSim" "not built"; return; }

    # --- fixed length ------------------------------------------------
    run frag-fixed "$prog" -n 200 -l 35 --seed 1 "$REF"
    ok "-n 200 -l 35 runs" "$RC" "$(lasterr frag-fixed)"
    cp "$W/frag-fixed.out" "$W/frag35.fa"
    eq "-n 200 produces 200 fragments" 200 "$(fa_count "$W/frag35.fa")"
    eq "-l 35 produces fragments of length 35 only" "35," "$(fa_lenset "$W/frag35.fa")"

    # the defline is chr:strand:start:end:length, the last field must agree
    # with the sequence that follows it
    local badlen
    badlen=$(awk '/^>/{n=split(substr($0,2),a,":"); want=a[n]; getline s;
		      if(length(s)!=want) bad++ } END{print bad+0}' "$W/frag35.fa")
    eq "the length in the defline matches the sequence" 0 "$badlen"

    local badcoord
    badcoord=$(awk '/^>/{n=split(substr($0,2),a,":");
			 if(a[n]!=(a[n-1]-a[n-2])) bad++ } END{print bad+0}' "$W/frag35.fa")
    eq "the coordinates in the defline span the fragment" 0 "$badcoord"

    eq "only ACGT is emitted" "" \
       "$(fa_seqs "$W/frag35.fa" | tr -d 'ACGT\n')"

    # --- the seed ----------------------------------------------------
    run frag-seedA "$prog" -n 50 -l 40 --seed 777 "$REF"
    run frag-seedB "$prog" -n 50 -l 40 --seed 777 "$REF"
    run frag-seedC "$prog" -n 50 -l 40 --seed 778 "$REF"
    if cmp -s "$W/frag-seedA.out" "$W/frag-seedB.out"; then
	pass "the same --seed gives the same fragments"
    else
	fail "the same --seed gives the same fragments"
    fi
    if cmp -s "$W/frag-seedA.out" "$W/frag-seedC.out"; then
	fail "a different --seed gives different fragments"
    else
	pass "a different --seed gives different fragments"
    fi

    # --- size distributions ------------------------------------------
    run frag-s "$prog" -n 300 -s "$W/sizes.size" --seed 2 "$REF"
    ok "-s reads a size distribution" "$RC" "$(lasterr frag-s)"
    eq "-s only draws lengths listed in the file" "30,40,50," \
       "$(fa_lenset "$W/frag-s.out")"

    run frag-f "$prog" -n 300 -f "$W/sizes.freq" --seed 3 "$REF"
    ok "-f reads a size frequency file" "$RC" "$(lasterr frag-f)"
    eq "-f only draws lengths listed in the file" "25,45," \
       "$(fa_lenset "$W/frag-f.out")"

    run frag-logn "$prog" -n 200 --loc 3.5 --scale 0.2 -m 20 -M 80 --seed 4 "$REF"
    ok "--loc/--scale draw from a lognormal" "$RC" "$(lasterr frag-logn)"
    local outofrange
    outofrange=$(fa_lens "$W/frag-logn.out" | awk '$1<20 || $1>80' | wc -l)
    eq "the lognormal stays within -m and -M" 0 "$outofrange"
    ne "the lognormal does not give a single length" 1 \
       "$(fa_lens "$W/frag-logn.out" | sort -u | wc -l)"

    run frag-mM "$prog" -n 300 -s "$W/sizes.size" -m 35 -M 45 --seed 5 "$REF"
    ok "-m/-M run with -s" "$RC" "$(lasterr frag-mM)"
    eq "-m 35 -M 45 keeps only the 40bp fragments" "40," "$(fa_lenset "$W/frag-mM.out")"

    # --- output formats ----------------------------------------------
    run frag-o "$prog" -n 100 -l 30 --seed 6 -o "$W/frag.fa.gz" "$REF"
    ok "-o writes a zipped fasta" "$RC" "$(lasterr frag-o)"
    if is_gzip "$W/frag.fa.gz"; then pass "-o output is gzipped"
    else fail "-o output is gzipped"; fi
    eq "-o keeps the requested number of fragments" 100 "$(fa_count "$W/frag.fa.gz")"

    run frag-b "$prog" -n 100 -l 30 --seed 6 -b "$W/frag.bam" "$REF"
    ok "-b writes a BAM" "$RC" "$(lasterr frag-b)"
    if is_bam "$W/frag.bam"; then pass "-b output is a BAM file"
    else fail "-b output is a BAM file" "no BAM magic"; fi
    if command -v samtools >/dev/null 2>&1; then
	eq "the BAM holds the requested number of fragments" 100 \
	   "$(samtools view -c "$W/frag.bam" 2>/dev/null)"
    else
	skip "the BAM holds the requested number of fragments" "no samtools"
    fi

    # -o and -b at the same time makes no sense
    run frag-ob "$prog" -n 10 -l 30 -o "$W/x.fa.gz" -b "$W/x.bam" "$REF"
    ne "-o together with -b is refused" 0 "$RC"

    # --- strand and naming -------------------------------------------
    run frag-norev "$prog" -n 100 -l 30 --norev --seed 8 "$REF"
    ok "--norev runs" "$RC" "$(lasterr frag-norev)"
    eq "--norev only reports the + strand" "+" \
       "$(grep '^>' "$W/frag-norev.out" | awk -F: '{print $2}' | sort -u | tr -d '\n')"
    eq "without --norev both strands are reported" "+-" \
       "$(grep '^>' "$W/frag35.fa" | awk -F: '{print $2}' | sort -u | tr -d '\n')"

    run frag-tag "$prog" -n 20 -l 30 --seed 9 -tag MYTAG "$REF"
    ok "-tag runs" "$RC" "$(lasterr frag-tag)"
    eq "-tag appends the tag to every defline" 20 \
       "$(grep -c '^>.*MYTAG$' "$W/frag-tag.out")"

    # -uniq used to swallow the argument that followed it
    run frag-uniq "$prog" -n 200 -l 30 --seed 10 -uniq "$REF"
    ok "-uniq runs and does not eat the next argument" "$RC" "$(lasterr frag-uniq)"
    eq "-uniq produces 200 fragments" 200 "$(fa_count "$W/frag-uniq.out")"
    eq "-uniq makes every defline unique" \
       "$(grep -c '^>' "$W/frag-uniq.out")" \
       "$(grep '^>' "$W/frag-uniq.out" | sort -u | wc -l | tr -d ' ')"

    # --- circular genomes --------------------------------------------
    run frag-circ "$prog" -n 200 -l 30 --seed 11 --circ testchr "$REF"
    ok "--circ runs" "$RC" "$(lasterr frag-circ)"
    eq "--circ produces 200 fragments" 200 "$(fa_count "$W/frag-circ.out")"

    # --- base composition --------------------------------------------
    # --comp takes the dnacomp.txt of mapDamage, not misincorporation.txt
    if [ -f "$DNACOMP" ]; then
	run frag-comp "$prog" -n 200 -l 30 --seed 12 --comp "$DNACOMP" --dist 5 "$REF"
	ok "--comp reads a mapDamage dnacomp file" "$RC" "$(lasterr frag-comp)"
	eq "--comp produces 200 fragments" 200 "$(fa_count "$W/frag-comp.out")"
	eq "--comp produces fragments of the requested length" "30," \
	   "$(fa_lenset "$W/frag-comp.out")"
    else
	skip "--comp reads a mapDamage dnacomp file" "no example file"
    fi

    # --- GC bias ------------------------------------------------------
    run frag-gc "$prog" -n 200 -l 30 --seed 13 -gc 0.5 "$REF"
    ok "-gc runs" "$RC" "$(lasterr frag-gc)"
    eq "-gc produces 200 fragments" 200 "$(fa_count "$W/frag-gc.out")"

    # --- fastq mode ---------------------------------------------------
    # --fq trims pre-selected fragments rather than drawing new ones
    awk 'BEGIN{
	     x=99;
	     for(r=1;r<=20;r++){
		 s="";
		 for(i=0;i<60;i++){ x=(16807*x)%2147483647; s=s substr("ACGT",int(x/1024)%4+1,1) }
		 printf "@read%d\n%s\n+\n%s\n", r, s, sprintf("%*s",60,"") ;
	     }
	 }' | sed 's/^ *$/IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII/' \
	    | gzip -c > "$W/reads.fq.gz"

    run frag-fq "$prog" --fq -l 40 --seed 14 "$W/reads.fq.gz"
    ok "--fq reads a fastq" "$RC" "$(lasterr frag-fq)"
    eq "--fq keeps every read" 20 "$(grep -c '^@read' "$W/frag-fq.out")"
    eq "--fq trims to the requested length" "40," \
       "$(awk 'NR%4==2{print length($0)}' "$W/frag-fq.out" | sort -n -u | tr '\n' ',')"
    eq "--fq keeps sequence and quality the same length" 0 \
       "$(awk 'NR%4==2{s=length($0)} NR%4==0{if(length($0)!=s) bad++} END{print bad+0}' "$W/frag-fq.out")"

    # --trim5p drops bases off the 5' end of the read and then takes the
    # requested length from there, so the read stays -l long
    run frag-trim "$prog" --fq -l 40 --trim5p 5 --seed 14 "$W/reads.fq.gz"
    ok "--trim5p runs" "$RC" "$(lasterr frag-trim)"
    eq "--trim5p keeps the reads at the requested length" "40," \
       "$(awk 'NR%4==2{print length($0)}' "$W/frag-trim.out" | sort -n -u | tr '\n' ',')"
    eq "--trim5p 5 drops the first 5 bases of the read" 0 \
       "$(paste <(gzip -cd "$W/reads.fq.gz" | awk 'NR%4==2') \
		<(awk 'NR%4==2' "$W/frag-trim.out") \
	  | awk '{ if(substr($1,6,40)!=$2) bad++ } END{print bad+0}')"
}

#####################################################################
# deamSim: adds deamination to fragments
#####################################################################

# counts the substitutions between two fastas holding the same records in the
# same order, prints "sub count" per line, e.g. "C>T 42"
subst_profile() {
    paste <(fa_seqs "$1") <(fa_seqs "$2") | awk '{
	n=length($1); if(length($2)<n) n=length($2);
	for(i=1;i<=n;i++){
	    a=toupper(substr($1,i,1)); b=toupper(substr($2,i,1));
	    if(a!=b) c[a">"b]++;
	}
    } END{ for(k in c) print k, c[k] }' | sort
}

# the fraction of bases that differ between two fastas holding the same records,
# counting only the positions [from..to] of each record, 1 based
diff_rate_range() {
    paste <(fa_seqs "$1") <(fa_seqs "$2") | awk -v f="$3" -v t="$4" '{
	n=length($1); e=(t<n)?t:n;
	for(i=f;i<=e;i++){
	    tot++;
	    if(toupper(substr($1,i,1))!=toupper(substr($2,i,1))) d++;
	}
    } END{ printf "%.6f\n", (tot?d/tot:0) }'
}

# The groups downstream of fragSim work on its output. Make that output on
# demand so that --only deamSim is just as complete as a full run.
ensure_fragments() {
    [ -x "$SRCDIR/fragSim" ] || return 1
    [ -s "$W/frag35.fa" ] || \
	"$SRCDIR/fragSim" -n 200 -l 35 --seed 1 "$REF" > "$W/frag35.fa" 2>/dev/null
    [ -s "$W/frag.bam" ] || \
	"$SRCDIR/fragSim" -n 100 -l 30 --seed 6 -b "$W/frag.bam" "$REF" >/dev/null 2>&1
    return 0
}

# greater [description] [a] [b]  -- passes when a > b
greater() {
    if awk -v a="$2" -v b="$3" 'BEGIN{exit !(a>b)}'; then pass "$1"
    else fail "$1" "expected $2 to be greater than $3"; fi
}

test_deamSim() {
    group "deamSim"
    local prog="$SRCDIR/deamSim"
    [ -x "$prog" ] || { skip "deamSim" "not built"; return; }
    ensure_fragments

    local nin
    nin=$(fa_count "$W/frag35.fa")

    # --- no damage at all --------------------------------------------
    run deam-none "$prog" -damage 0,0,0,0 --seed 21 "$W/frag35.fa"
    ok "-damage 0,0,0,0 runs" "$RC" "$(lasterr deam-none)"
    eq "-damage 0,0,0,0 leaves the fragments untouched" "" \
       "$(subst_profile "$W/frag35.fa" "$W/deam-none.out")"

    # --- the Briggs model ---------------------------------------------
    run deam-briggs "$prog" -damage 0.03,0.4,0.01,0.7 --seed 22 "$W/frag35.fa"
    ok "-damage runs the Briggs model" "$RC" "$(lasterr deam-briggs)"
    eq "-damage keeps every fragment" "$nin" "$(fa_count "$W/deam-briggs.out")"
    eq "-damage does not change any fragment length" "35," "$(fa_lenset "$W/deam-briggs.out")"
    eq "-damage only deaminates, i.e. C>T and G>A" "" \
       "$(subst_profile "$W/frag35.fa" "$W/deam-briggs.out" | awk '$1!="C>T" && $1!="G>A"')"
    ne "-damage actually introduces substitutions" "" \
       "$(subst_profile "$W/frag35.fa" "$W/deam-briggs.out")"

    run deam-briggs2 "$prog" -damage 0.03,0.4,0.01,0.7 --seed 22 "$W/frag35.fa"
    if cmp -s "$W/deam-briggs.out" "$W/deam-briggs2.out"; then
	pass "the same --seed gives the same deamination"
    else
	fail "the same --seed gives the same deamination"
    fi
    run deam-briggs3 "$prog" -damage 0.03,0.4,0.01,0.7 --seed 23 "$W/frag35.fa"
    if cmp -s "$W/deam-briggs.out" "$W/deam-briggs3.out"; then
	fail "a different --seed gives a different deamination"
    else
	pass "a different --seed gives a different deamination"
    fi

    # --- the built in matrices ----------------------------------------
    run deam-single "$prog" -mat single --seed 24 "$W/frag35.fa"
    ok "-mat single runs" "$RC" "$(lasterr deam-single)"
    eq "-mat single keeps every fragment" "$nin" "$(fa_count "$W/deam-single.out")"

    run deam-double "$prog" -mat double --seed 25 "$W/frag35.fa"
    ok "-mat double runs" "$RC" "$(lasterr deam-double)"
    eq "-mat double keeps every fragment" "$nin" "$(fa_count "$W/deam-double.out")"

    eq "-mat only deaminates, i.e. C>T and G>A" "" \
       "$(subst_profile "$W/frag35.fa" "$W/deam-double.out" | awk '$1!="C>T" && $1!="G>A"')"

    # the point of the whole thing: damage must pile up at the ends of the
    # molecule and thin out in the middle.  The built in matrices are empirical
    # so both C>T and G>A occur at both ends, only the rate carries the signal,
    # which needs more than a couple of hundred fragments to come out cleanly.
    "$SRCDIR/fragSim" -n 5000 -l 60 --seed 33 "$REF" > "$W/frag60.fa" 2>/dev/null
    local mid end3 end5 mat
    for mat in single double; do
	run "deam-rate-$mat" "$prog" -mat "$mat" --seed 34 "$W/frag60.fa"
	ok "-mat $mat runs on 5000 fragments" "$RC" "$(lasterr "deam-rate-$mat")"
	end5=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-$mat.out" 1 3)
	end3=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-$mat.out" 58 60)
	mid=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-$mat.out" 25 36)
	greater "-mat $mat damages the 5' end more than the middle" "$end5" "$mid"
	greater "-mat $mat damages the 3' end more than the middle" "$end3" "$mid"
    done

    run deam-rate-briggs "$prog" -damage 0.03,0.4,0.01,0.7 --seed 35 "$W/frag60.fa"
    ok "-damage runs on 5000 fragments" "$RC" "$(lasterr deam-rate-briggs)"
    end5=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-briggs.out" 1 3)
    end3=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-briggs.out" 58 60)
    mid=$(diff_rate_range "$W/frag60.fa" "$W/deam-rate-briggs.out" 25 36)
    greater "-damage puts more damage on the 5' end than in the middle" "$end5" "$mid"
    greater "-damage puts more damage on the 3' end than in the middle" "$end3" "$mid"

    # a stronger single strand deamination probability must give more damage
    run deam-weak   "$prog" -damage 0.03,0.4,0.01,0.1 --seed 36 "$W/frag60.fa"
    run deam-strong "$prog" -damage 0.03,0.4,0.01,0.9 --seed 36 "$W/frag60.fa"
    greater "a higher single strand deamination rate gives more substitutions" \
	    "$(diff_rate_range "$W/frag60.fa" "$W/deam-strong.out" 1 60)" \
	    "$(diff_rate_range "$W/frag60.fa" "$W/deam-weak.out" 1 60)"

    # --- matrix files -------------------------------------------------
    if [ -f "$SRCDIR/matrices/single-5.dat" ]; then
	run deam-matfile "$prog" -matfile "$SRCDIR/matrices/single-" --seed 26 "$W/frag35.fa"
	ok "-matfile reads src/matrices" "$RC" "$(lasterr deam-matfile)"
	eq "-matfile keeps every fragment" "$nin" "$(fa_count "$W/deam-matfile.out")"
    else
	skip "-matfile reads src/matrices" "no matrices shipped"
    fi

    # --- a mapDamage misincorporation file -----------------------------
    if [ -f "$MAPDAMAGE" ]; then
	run deam-md "$prog" -mapdamage "$MAPDAMAGE" double --seed 27 "$W/frag35.fa"
	ok "-mapdamage double runs" "$RC" "$(lasterr deam-md)"
	eq "-mapdamage keeps every fragment" "$nin" "$(fa_count "$W/deam-md.out")"
	eq "-mapdamage only deaminates" "" \
	   "$(subst_profile "$W/frag35.fa" "$W/deam-md.out" | awk '$1!="C>T" && $1!="G>A"')"

	run deam-mds "$prog" -mapdamage "$MAPDAMAGE" single --seed 28 "$W/frag35.fa"
	ok "-mapdamage single runs" "$RC" "$(lasterr deam-mds)"
    else
	skip "-mapdamage" "no example file"
    fi

    # --- a .prof file, as produced by mapDamage2prof --------------------
    if [ -x "$SRCDIR/mapDamage2prof" ] && [ -f "$MAPDAMAGE" ]; then
	"$SRCDIR/mapDamage2prof" -double -5p "$W/p5p.prof" -3p "$W/p3p.prof" "$MAPDAMAGE" \
	    >/dev/null 2>&1
	run deam-prof "$prog" -profile "$W/p" --seed 29 "$W/frag35.fa"
	ok "-profile reads the output of mapDamage2prof" "$RC" "$(lasterr deam-prof)"
	eq "-profile keeps every fragment" "$nin" "$(fa_count "$W/deam-prof.out")"
	eq "-profile only deaminates" "" \
	   "$(subst_profile "$W/frag35.fa" "$W/deam-prof.out" | awk '$1!="C>T" && $1!="G>A"')"
    else
	skip "-profile" "mapDamage2prof or the example file is missing"
    fi

    # --- output formats -------------------------------------------------
    run deam-o "$prog" -damage 0.03,0.4,0.01,0.7 --seed 30 -o "$W/deam.fa.gz" "$W/frag35.fa"
    ok "-o writes a zipped fasta" "$RC" "$(lasterr deam-o)"
    if is_gzip "$W/deam.fa.gz"; then pass "-o output is gzipped"
    else fail "-o output is gzipped"; fi
    eq "-o keeps every fragment" "$nin" "$(fa_count "$W/deam.fa.gz")"

    if [ -s "$W/frag.bam" ]; then
	run deam-bam "$prog" -damage 0.03,0.4,0.01,0.7 --seed 31 -b "$W/deam.bam" "$W/frag.bam"
	ok "-b reads and writes a BAM" "$RC" "$(lasterr deam-bam)"
	if is_bam "$W/deam.bam"; then pass "-b output is a BAM file"
	else fail "-b output is a BAM file" "no BAM magic"; fi
    else
	skip "-b reads and writes a BAM" "fragSim produced no BAM"
    fi

    run deam-name "$prog" -damage 0.03,0.4,1.0,1.0 --seed 32 -name "$W/frag35.fa"
    ok "-name runs" "$RC" "$(lasterr deam-name)"
    ne "-name records the deaminated bases in the defline" 0 \
       "$(grep -c '^>.*DEAM' "$W/deam-name.out")"
}

#####################################################################
# adptSim: turns fragments into reads with adapters
#####################################################################

revcomp() {
    awk '{ s=""; for(i=length($0);i>0;i--) s=s substr($0,i,1); print s }' \
	| tr 'ACGTacgt' 'TGCAtgca'
}

test_adptSim() {
    group "adptSim"
    local prog="$SRCDIR/adptSim"
    [ -x "$prog" ] || { skip "adptSim" "not built"; return; }
    ensure_fragments

    local nin
    nin=$(fa_count "$W/frag35.fa")

    # --- paired end, zipped fasta ---------------------------------------
    run adpt-fr "$prog" -fr "$W/fwd.fa.gz" -rr "$W/rev.fa.gz" -l 50 "$W/frag35.fa"
    ok "-fr/-rr run" "$RC" "$(lasterr adpt-fr)"
    eq "-fr keeps every fragment" "$nin" "$(fa_count "$W/fwd.fa.gz")"
    eq "-rr keeps every fragment" "$nin" "$(fa_count "$W/rev.fa.gz")"
    eq "-l 50 gives forward reads of length 50" "50," "$(fa_lenset "$W/fwd.fa.gz")"
    eq "-l 50 gives reverse reads of length 50" "50," "$(fa_lenset "$W/rev.fa.gz")"

    # the read must start with the fragment and only then run into the adapter
    eq "the forward read starts with the fragment" 0 \
       "$(paste <(fa_seqs "$W/frag35.fa") <(fa_seqs "$W/fwd.fa.gz") \
	  | awk '{ if(substr($2,1,length($1))!=toupper($1)) bad++ } END{print bad+0}')"
    eq "the reverse read starts with the reverse complement of the fragment" 0 \
       "$(paste <(fa_seqs "$W/frag35.fa" | revcomp) <(fa_seqs "$W/rev.fa.gz") \
	  | awk '{ if(substr($2,1,length($1))!=toupper($1)) bad++ } END{print bad+0}')"

    # --- custom adapters -------------------------------------------------
    local ADPF="TTTTTTTTTTTTTTTTTTTT" ADPR="GGGGGGGGGGGGGGGGGGGG"
    run adpt-custom "$prog" -f "$ADPF" -s "$ADPR" \
	-fr "$W/fwd2.fa.gz" -rr "$W/rev2.fa.gz" -l 50 "$W/frag35.fa"
    ok "-f/-s set the adapters" "$RC" "$(lasterr adpt-custom)"
    eq "the forward adapter follows the fragment" 0 \
       "$(fa_seqs "$W/fwd2.fa.gz" | awk '{ if(substr($0,36)!="TTTTTTTTTTTTTTT") bad++ } END{print bad+0}')"
    eq "the reverse adapter follows the fragment" 0 \
       "$(fa_seqs "$W/rev2.fa.gz" | awk '{ if(substr($0,36)!="GGGGGGGGGGGGGGG") bad++ } END{print bad+0}')"

    # --- a read shorter than the fragment: no adapter, just a truncation ---
    run adpt-short "$prog" -fr "$W/fwd3.fa.gz" -rr "$W/rev3.fa.gz" -l 20 "$W/frag35.fa"
    ok "a read length below the fragment length runs" "$RC" "$(lasterr adpt-short)"
    eq "the fragment is truncated to the read length" "20," "$(fa_lenset "$W/fwd3.fa.gz")"
    eq "no adapter is added when the fragment is long enough" 0 \
       "$(paste <(fa_seqs "$W/frag35.fa") <(fa_seqs "$W/fwd3.fa.gz") \
	  | awk '{ if(substr(toupper($1),1,20)!=$2) bad++ } END{print bad+0}')"

    # --- the ART flavours -------------------------------------------------
    run adpt-arts "$prog" -arts "$W/art_s.fa" -l 50 "$W/frag35.fa"
    ok "-arts writes single end ART input" "$RC" "$(lasterr adpt-arts)"
    eq "-arts keeps every fragment" "$nin" "$(fa_count "$W/art_s.fa")"
    eq "-arts gives reads of length 50" "50," "$(fa_lenset "$W/art_s.fa")"

    run adpt-artp "$prog" -artp "$W/art_p.fa" -l 50 "$W/frag35.fa"
    ok "-artp writes wrapped around paired end ART input" "$RC" "$(lasterr adpt-artp)"
    eq "-artp keeps every fragment" "$nin" "$(fa_count "$W/art_p.fa")"
    eq "-artp wraps the two mates into one record of twice the read length" "100," \
       "$(fa_lenset "$W/art_p.fa")"

    # --- BAM --------------------------------------------------------------
    if [ -s "$W/frag.bam" ]; then
	run adpt-bs "$prog" -bs "$W/adpt_s.bam" -l 50 "$W/frag.bam"
	ok "-bs writes a single end BAM" "$RC" "$(lasterr adpt-bs)"
	if is_bam "$W/adpt_s.bam"; then pass "-bs output is a BAM file"
	else fail "-bs output is a BAM file" "no BAM magic"; fi

	run adpt-bp "$prog" -bp "$W/adpt_p.bam" -l 50 "$W/frag.bam"
	ok "-bp writes a paired end BAM" "$RC" "$(lasterr adpt-bp)"
	if is_bam "$W/adpt_p.bam"; then pass "-bp output is a BAM file"
	else fail "-bp output is a BAM file" "no BAM magic"; fi

	if command -v samtools >/dev/null 2>&1; then
	    eq "-bp writes two records per fragment" \
	       "$(( $(samtools view -c "$W/adpt_s.bam" 2>/dev/null) * 2 ))" \
	       "$(samtools view -c "$W/adpt_p.bam" 2>/dev/null)"
	else
	    skip "-bp writes two records per fragment" "no samtools"
	fi
    else
	skip "-bs/-bp" "fragSim produced no BAM"
    fi

    # --- tagging -----------------------------------------------------------
    run adpt-tag "$prog" -tag _MYTAG -fr "$W/fwd4.fa.gz" -rr "$W/rev4.fa.gz" -l 50 "$W/frag35.fa"
    ok "-tag runs" "$RC" "$(lasterr adpt-tag)"
    eq "-tag appends the tag to every defline" "$nin" \
       "$(fa_cat "$W/fwd4.fa.gz" | grep -c '^>.*_MYTAG$')"
}

#####################################################################
# mapDamage2prof: mapDamage misincorporation file -> .prof
#####################################################################

# checks that a .prof file looks the way deamSim expects it to
check_prof() {
    local file="$1" label="$2"
    if [ ! -s "$file" ]; then fail "$label is written"; return 1; fi
    pass "$label is written"

    eq "$label has the 12 substitution header" \
       "A>C	A>G	A>T	C>A	C>G	C>T	G>A	G>C	G>T	T>A	T>C	T>G" \
       "$(head -1 "$file")"

    ne "$label has at least one position" 0 "$(( $(wc -l < "$file") - 1 ))"

    eq "$label has 12 columns on every row" 0 \
       "$(tail -n +2 "$file" | awk -F'\t' 'NF!=12{bad++} END{print bad+0}')"

    eq "$label holds probabilities between 0 and 1" 0 \
       "$(tail -n +2 "$file" | awk -F'\t' '{
	      for(i=1;i<=NF;i++){ if($i+0<0 || $i+0>1) bad++ }
	  } END{print bad+0}')"
    return 0
}

# is column [n] of the first position of a .prof clearly above zero?
prof_nonzero() {
    awk -F'\t' -v c="$2" 'NR==2{print ($c+0 > 0.001)?"yes":"no"}' "$1"
}

# how many cells outside the C>T and G>A columns are not zero? deamSim is only
# ever handed deamination, everything else in a .prof has to be zeroed out
prof_other_than_deamination() {
    awk -F'\t' 'FNR>1{
	for(i=1;i<=NF;i++){ if(i!=6 && i!=7 && $i+0!=0) bad++ }
    } END{print bad+0}' "$@"
}

test_mapDamage2prof() {
    group "mapDamage2prof"
    local prog="$SRCDIR/mapDamage2prof"
    [ -x "$prog" ] || { skip "mapDamage2prof" "not built"; return; }
    [ -f "$MAPDAMAGE" ] || { skip "mapDamage2prof" "no example file"; return; }

    run md-double "$prog" -double -5p "$W/md5p.prof" -3p "$W/md3p.prof" "$MAPDAMAGE"
    ok "-double runs on the LaBrana example" "$RC" "$(lasterr md-double)"
    check_prof "$W/md5p.prof" "the -double 5' profile"
    check_prof "$W/md3p.prof" "the -double 3' profile"

    # column 6 is C>T and column 7 is G>A.  The 5' end always reports C>T only,
    # what tells the three modes apart is the 3' end.
    eq "-double reports C>T at the first 5' position" "yes" "$(prof_nonzero "$W/md5p.prof" 6)"
    eq "-double leaves G>A out of the 5' profile" "no"  "$(prof_nonzero "$W/md5p.prof" 7)"
    eq "-double reports G>A at the first 3' position" "yes" "$(prof_nonzero "$W/md3p.prof" 7)"
    eq "-double leaves C>T out of the 3' profile" "no"  "$(prof_nonzero "$W/md3p.prof" 6)"
    eq "-double zeroes every substitution that is not deamination" 0 \
       "$(prof_other_than_deamination "$W/md5p.prof" "$W/md3p.prof")"

    run md-single "$prog" -single -5p "$W/ms5p.prof" -3p "$W/ms3p.prof" "$MAPDAMAGE"
    ok "-single runs" "$RC" "$(lasterr md-single)"
    check_prof "$W/ms5p.prof" "the -single 5' profile"
    check_prof "$W/ms3p.prof" "the -single 3' profile"
    # a single stranded library carries C>T at both ends
    eq "-single reports C>T at the first 3' position" "yes" "$(prof_nonzero "$W/ms3p.prof" 6)"
    eq "-single leaves G>A out of the 3' profile" "no"  "$(prof_nonzero "$W/ms3p.prof" 7)"

    run md-both "$prog" -both -5p "$W/mb5p.prof" -3p "$W/mb3p.prof" "$MAPDAMAGE"
    ok "-both runs" "$RC" "$(lasterr md-both)"
    check_prof "$W/mb5p.prof" "the -both 5' profile"
    check_prof "$W/mb3p.prof" "the -both 3' profile"
    eq "-both reports C>T at the 3' end" "yes" "$(prof_nonzero "$W/mb3p.prof" 6)"
    eq "-both reports G>A at the 3' end too" "yes" "$(prof_nonzero "$W/mb3p.prof" 7)"

    # the deamination rate has to decay away from the end of the molecule
    greater "the C>T rate decays along the 5' end" \
	    "$(awk -F'\t' 'NR==2{print $6}' "$W/md5p.prof")" \
	    "$(awk -F'\t' 'NR==6{print $6}' "$W/md5p.prof")"
    greater "the G>A rate decays along the 3' end" \
	    "$(awk -F'\t' 'NR==2{print $7}' "$W/md3p.prof")" \
	    "$(awk -F'\t' 'NR==6{print $7}' "$W/md3p.prof")"

    run md-h "$prog" -h -double -5p "$W/mh5p.prof" -3p "$W/mh3p.prof" "$MAPDAMAGE"
    ok "-h runs" "$RC" "$(lasterr md-h)"
    eq "-h adds a pos column to the header" "pos" \
       "$(head -1 "$W/mh5p.prof" | cut -f1)"

    if [ -f "$MAPDAMAGE2" ]; then
	run md-ust "$prog" -double -5p "$W/mu5p.prof" -3p "$W/mu3p.prof" "$MAPDAMAGE2"
	ok "-double runs on the Ust Ishim example" "$RC" "$(lasterr md-ust)"
	check_prof "$W/mu5p.prof" "the Ust Ishim 5' profile"
    else
	skip "the Ust Ishim example" "no example file"
    fi

    run md-missing "$prog" -double -5p "$W/x.prof" -3p "$W/y.prof" "$W/does-not-exist.txt"
    ne "a missing input file is refused" 0 "$RC"
}

#####################################################################
# damage_patterns2prof: damage-patterns matrices -> .prof
#####################################################################

test_damage_patterns2prof() {
    group "damage_patterns2prof"
    local prog="$SRCDIR/damage_patterns2prof"
    [ -x "$prog" ] || { skip "damage_patterns2prof" "not built"; return; }

    run dp-double "$prog" -double -5p "$W/dp5p.prof" -3p "$W/dp3p.prof" \
	"$W/dp-5p.dat" "$W/dp-3p.dat"
    ok "-double runs" "$RC" "$(lasterr dp-double)"
    check_prof "$W/dp5p.prof" "the -double 5' profile"
    check_prof "$W/dp3p.prof" "the -double 3' profile"

    eq "the number of positions matches the input matrix" \
       "$(( $(wc -l < "$W/dp-5p.dat") - 1 ))" \
       "$(( $(wc -l < "$W/dp5p.prof") - 1 ))"

    # the fixture puts C>T at the 5' end and G>A at the 3' end
    eq "the 5' profile picks up C>T at the end of the molecule" "yes" \
       "$(prof_nonzero "$W/dp5p.prof" 6)"
    eq "the 3' profile picks up G>A at the end of the molecule" "yes" \
       "$(prof_nonzero "$W/dp3p.prof" 7)"

    # the 5' matrix is copied straight through, position for position
    eq "the 5' profile keeps the order of the input matrix" 0 \
       "$(paste <(tail -n +2 "$W/dp5p.prof" | cut -f6) \
		<(tail -n +2 "$W/dp-5p.dat" | cut -f7 | awk '{print $1}') \
	  | awk '{ d=$1-$2; if(d<0) d=-d; if(d>1e-6) bad++ } END{print bad+0}')"

    # damage-patterns writes the 3' matrix the other way round, so the tool has
    # to flip it before deamSim can read it
    eq "the 3' profile is flipped so row 1 is the base next to the 3' end" 0 \
       "$(paste <(tail -n +2 "$W/dp3p.prof" | cut -f7) \
		<(tail -n +2 "$W/dp-3p.dat" | cut -f8 | awk '{print $1}' | tac) \
	  | awk '{ d=$1-$2; if(d<0) d=-d; if(d>1e-6) bad++ } END{print bad+0}')"

    run dp-both "$prog" -both -5p "$W/dpb5p.prof" -3p "$W/dpb3p.prof" \
	"$W/dp-5p.dat" "$W/dp-3p.dat"
    ok "-both runs" "$RC" "$(lasterr dp-both)"
    check_prof "$W/dpb5p.prof" "the -both 5' profile"

    run dp-single "$prog" -single -5p "$W/dps5p.prof" -3p "$W/dps3p.prof" \
	"$W/dp-5p.dat" "$W/dp-3p.dat"
    ok "-single runs" "$RC" "$(lasterr dp-single)"
    check_prof "$W/dps5p.prof" "the -single 5' profile"

    # unlike mapDamage2prof, this tool passes the whole matrix through and
    # -single/-double/-both make no difference to what it writes
    if cmp -s "$W/dp5p.prof" "$W/dps5p.prof" && cmp -s "$W/dp5p.prof" "$W/dpb5p.prof" \
	    && cmp -s "$W/dp3p.prof" "$W/dps3p.prof" && cmp -s "$W/dp3p.prof" "$W/dpb3p.prof"; then
	pass "-single, -double and -both write the same profile"
    else
	fail "-single, -double and -both write the same profile"
    fi

    # a header that is not the one damage-patterns writes must be rejected
    printf 'not\ta\theader\n0\t1\t2\n' > "$W/bad.dat"
    run dp-bad "$prog" -5p "$W/x.prof" -3p "$W/y.prof" "$W/bad.dat" "$W/bad.dat"
    ne "a malformed matrix is refused" 0 "$RC"

    # the output must be readable by deamSim
    ensure_fragments
    if [ -x "$SRCDIR/deamSim" ] && [ -s "$W/frag35.fa" ]; then
	cp "$W/dp5p.prof" "$W/dpprof5p.prof"
	cp "$W/dp3p.prof" "$W/dpprof3p.prof"
	run dp-deam "$SRCDIR/deamSim" -profile "$W/dpprof" --seed 41 "$W/frag35.fa"
	ok "deamSim reads the resulting profile" "$RC" "$(lasterr dp-deam)"
	ne "the resulting profile deaminates something" "" \
	   "$(subst_profile "$W/frag35.fa" "$W/dp-deam.out")"
    else
	skip "deamSim reads the resulting profile" "deamSim is missing"
    fi
}

#####################################################################
# pipeline: the three tools chained the way gargammel.pl chains them
#####################################################################

test_pipeline() {
    group "pipeline"
    for p in fragSim deamSim adptSim; do
	[ -x "$SRCDIR/$p" ] || { skip "pipeline" "$p is not built"; return; }
    done

    "$SRCDIR/fragSim" -n 500 -s "$W/sizes.size" --seed 51 "$REF" > "$W/pipe1.fa" 2>"$W/pipe1.err"
    local rc1=$?
    "$SRCDIR/deamSim" -mat double --seed 52 "$W/pipe1.fa" > "$W/pipe2.fa" 2>"$W/pipe2.err"
    local rc2=$?
    "$SRCDIR/adptSim" -fr "$W/pipe_f.fa.gz" -rr "$W/pipe_r.fa.gz" -l 75 "$W/pipe2.fa" \
	>"$W/pipe3.out" 2>"$W/pipe3.err"
    local rc3=$?

    ok "fragSim | deamSim | adptSim runs end to end" "$(( rc1 + rc2 + rc3 ))"
    eq "no fragment is lost between fragSim and deamSim" \
       "$(fa_count "$W/pipe1.fa")" "$(fa_count "$W/pipe2.fa")"
    eq "no fragment is lost between deamSim and adptSim" \
       "$(fa_count "$W/pipe2.fa")" "$(fa_count "$W/pipe_f.fa.gz")"
    eq "every read comes out at the requested length" "75," "$(fa_lenset "$W/pipe_f.fa.gz")"
    eq "deamSim preserves the fragment lengths" \
       "$(fa_lenset "$W/pipe1.fa")" "$(fa_lenset "$W/pipe2.fa")"
    eq "the deflines survive the whole pipeline" \
       "$(grep -c '^>' "$W/pipe1.fa")" \
       "$(paste <(grep '^>' "$W/pipe1.fa") <(fa_cat "$W/pipe_f.fa.gz" | grep '^>') \
	  | awk '$1==$2' | wc -l | tr -d ' ')"
}

#####################################################################
# gargammel.pl: the wrapper
#####################################################################

test_gargammel() {
    group "gargammel.pl"
    local prog="$GARGDIR/gargammel.pl"
    [ -x "$prog" ] || { skip "gargammel.pl" "not found"; return; }
    if ! command -v perl >/dev/null 2>&1; then skip "gargammel.pl" "no perl"; return; fi

    run garg-help perl "$prog" --help
    ne "--help prints something" "" "$(cat "$W/garg-help.out" "$W/garg-help.err")"

    # the wrapper checks every executable it drives up front, so without ART
    # even --mock refuses to start
    if [ ! -x "$GARGDIR/art_src_MountRainier/art_illumina" ]; then
	skip "the gargammel.pl pipeline" "art_illumina is not built, run 'make'"
	return
    fi

    # the input directory the wrapper expects
    local IN="$W/inputfolder"
    mkdir -p "$IN/endo" "$IN/cont" "$IN/bact"
    make_fasta "$IN/endo/endo.1.fa" testchr 20000 42
    make_fasta "$IN/endo/endo.2.fa" testchr 20000 43
    make_fasta "$IN/cont/cont.1.fa" testchr 20000 44
    make_fasta "$IN/cont/cont.2.fa" testchr 20000 45
    make_fasta "$IN/bact/bact.1.fa" bactchr 20000 46
    printf 'bact.1.fa\t1.0\n' > "$IN/bact/list"

    run garg-mock perl "$prog" --mock --comp 0,0,1 -n 1000 -l 35 \
	-damage 0.03,0.4,0.01,0.3 -o "$W/mockout" "$IN"
    ok "--mock runs" "$RC" "$(lasterr garg-mock)"
    if grep -q 'fragSim' "$W/garg-mock.out" "$W/garg-mock.err"; then
	pass "--mock reports the fragSim command"
    else
	fail "--mock reports the fragSim command"
    fi
    if grep -q 'deamSim' "$W/garg-mock.out" "$W/garg-mock.err"; then
	pass "--mock reports the deamSim command"
    else
	fail "--mock reports the deamSim command"
    fi
    if grep -q 'adptSim' "$W/garg-mock.out" "$W/garg-mock.err"; then
	pass "--mock reports the adptSim command"
    else
	fail "--mock reports the adptSim command"
    fi
    if [ ! -s "$W/mockout_s1.fq.gz" ] && [ ! -s "$W/mockout.fq.gz" ]; then
	pass "--mock writes no reads"
    else
	fail "--mock writes no reads" "it produced output anyway"
    fi

    run garg-full perl "$prog" --comp 0,0.1,0.9 -n 500 -l 40 \
	-damage 0.03,0.4,0.01,0.3 -rl 60 -o "$W/sim" "$IN"
    ok "a full paired end simulation runs" "$RC" "$(lasterr garg-full)"
    if [ -s "$W/sim_s1.fq.gz" ] && [ -s "$W/sim_s2.fq.gz" ]; then
	pass "both mates are written"
	eq "the two mates hold the same number of reads" \
	   "$(gzip -cd "$W/sim_s1.fq.gz" | wc -l)" \
	   "$(gzip -cd "$W/sim_s2.fq.gz" | wc -l)"
	eq "the reads come out at the requested length" "60," \
	   "$(gzip -cd "$W/sim_s1.fq.gz" | awk 'NR%4==2{print length($0)}' | sort -n -u | tr '\n' ',')"
	eq "sequence and quality have the same length" 0 \
	   "$(gzip -cd "$W/sim_s1.fq.gz" \
	      | awk 'NR%4==2{s=length($0)} NR%4==0{if(length($0)!=s) bad++} END{print bad+0}')"
    else
	fail "both mates are written" "no fastq in $W"
    fi

    run garg-se perl "$prog" --comp 0,0,1 -n 200 -l 40 -se -rl 60 -o "$W/simse" "$IN"
    ok "a single end simulation runs" "$RC" "$(lasterr garg-se)"
    if [ -s "$W/simse_s.fq.gz" ]; then
	pass "single end reads are written"
	eq "-se writes one file, not a pair" "" \
	   "$(ls "$W"/simse_s1.fq.gz "$W"/simse_s2.fq.gz 2>/dev/null)"
	eq "-se produces one read per fragment" 200 \
	   "$(( $(gzip -cd "$W/simse_s.fq.gz" | wc -l) / 4 ))"
	eq "-se honours the read length" "60," \
	   "$(gzip -cd "$W/simse_s.fq.gz" | awk 'NR%4==2{print length($0)}' | sort -n -u | tr '\n' ',')"
    else
	fail "single end reads are written" "no fastq in $W"
    fi
}

#####################################################################

printf '%sgargammel test suite%s\n' "$CBOLD" "$COFF"
printf 'repository:  %s\n' "$GARGDIR"
printf 'working dir: %s%s\n' "$W" "$([ -n "$KEEP" ] && echo ' (kept)')"

setup_fixtures

for g in $TESTGROUPS; do
    wants "$g" || continue
    "test_$g"
done

printf '\n%s== summary%s\n' "$CBOLD" "$COFF"
printf '  %d passed, %d failed, %d skipped\n' "$NPASS" "$NFAIL" "$NSKIP"

if [ "$NFAIL" -gt 0 ]; then
    printf '\n%sfailed tests:%s%s\n' "$CRED" "$COFF" "$FAILED"
    exit 1
fi
exit 0
