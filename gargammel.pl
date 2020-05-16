#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use File::Copy;
use Time::HiRes qw/ time sleep /;

my $mock =0;


sub checkDeamParam{
  my ($briggs,$matfile) = @_;


  if ( (defined $matfile) ) {

    if ( !(-e $matfile."5.dat")) {
      die "Matrix file ".$matfile."5.dat does not exists\n";
    }

    if ( !(-e $matfile."3.dat")) {
      die "Matrix file ".$matfile."3.dat does not exists\n";
    }

  }

  if ( (defined $briggs) ) {
    my @arr=split(",",$briggs);

    if ($#arr != 3) {
      die "The option for -damage must be 4 comma-delimited numbers, found ".($#arr+1)." parameters in ".$briggs."\n";
    }

    if (!looks_like_number($arr[0])) {
      die "The option for -damage must be 4 comma-delimited numbers, parameter 1 :".$arr[0]." must be a number\n";
    }
    if (!looks_like_number($arr[1])) {
      die "The option for -damage must be 4 comma-delimited numbers, parameter 2 :".$arr[1]." must be a number\n";
    }
    if (!looks_like_number($arr[2])) {
      die "The option for -damage must be 4 comma-delimited numbers, parameter 3 :".$arr[2]." must be a number\n";
    }
    if (!looks_like_number($arr[3])) {
      die "The option for -damage must be 4 comma-delimited numbers, parameter 4 :".$arr[3]." must be a number\n";
    }
  }


  if ( (defined $matfile) &&
       (defined $briggs)  ) {
    die "Specify either -matfile or -damage but not both\n";
  }

}

sub copycmd{
  my ($source,$destination) = @_;
  print STDERR "copying  ". $source." to ".$destination."\n";
  if($mock != 1){
    copy( $source,$destination ) or die "Copy file failed: $!";
  }
}


sub runcmd{
  my ($cmdtorun) = @_;

  print STDERR "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    my @argstorun = ( "bash", "-c", $cmdtorun );

    if(system(@argstorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}

sub runcmdforce{
  my ($cmdtorun) = @_;

  print STDERR "running cmd ". $cmdtorun."\n";

  my @argstorun = ( "bash", "-c", $cmdtorun );

  if (system(@argstorun) != 0) {
    die "system  cmd $cmdtorun failed: $?"
  } else {
  }
}






sub runcmdReturnOutput {
  my $command = join ' ', @_;
  reverse ($_ = qx{$command 2>&1}, $? >> 8);
}

sub listdir {
  my ($dirtolist) = @_;
  my @arrayToReturn;
  opendir(DIR,$dirtolist) or die "Cannot open directory ".$!;
  while(my $dir = readdir(DIR)){
    if($dir eq "."){
      next;
    }

    if($dir eq ".."){
      next;
    }

    if(-d $dirtolist."".$dir){
      #print "dir ".$dirtolist."".$dir."\n";
      push(@arrayToReturn,$dirtolist."".$dir."/");
    }
  }
  close(DIR);
  return @arrayToReturn;
}


sub listdirFai {
  my ($dirtolist) = @_;
  my @arrayToReturn;
  opendir(DIR,$dirtolist) or die "Cannot open directory ".$!;
  while(my $dir = readdir(DIR)){
    if($dir eq "."){
      next;
    }

    if($dir eq ".."){
      next;
    }

    if(-f $dirtolist."".$dir){
      #print "dir ".$dirtolist."".$dir."\n";
      if($dir =~ /.fai$/  ){
	push(@arrayToReturn,$dirtolist."".$dir);
      }
    }
  }
  close(DIR);
  return @arrayToReturn;
}

sub listdirFa {
  my ($dirtolist) = @_;
  my @arrayToReturn;

  opendir(DIR,$dirtolist) or die "Cannot open directory ".$!;
  while(my $dir = readdir(DIR)){
    if($dir eq "."){
      next;
    }

    if($dir eq ".."){
      next;
    }

    if(-f $dirtolist."".$dir){

      if($dir =~ /.fa$/       ||
	 $dir =~ /.fa.gz$/    ||
	 $dir =~ /.fna$/      ||
	 $dir =~ /.fasta$/    ||
	 $dir =~ /.fasta.gz$/ ){
	push(@arrayToReturn,$dirtolist."".$dir);
      }
    }
  }
  close(DIR);
  return @arrayToReturn;
}


sub listdirCdir {
  my ($dirtolist) = @_;
  my @arrayToReturn;
  opendir(DIR,$dirtolist) or die "Cannot open directory ".$!;
  while(my $dir = readdir(DIR)){
    if($dir eq "."){
      next;
    }

    if($dir eq ".."){
      next;
    }

    if(-d $dirtolist."".$dir){
      #print "dir ".$dirtolist."".$dir."\n";
      #if($dir =~ /^C[0-9]+$/        ){
      push(@arrayToReturn,$dirtolist."".$dir);
      #}else{


    }
  }
  close(DIR);
  return @arrayToReturn;
}


sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}


sub isZipped{
  my ($fileToTest) = @_;
  if(runcmdReturnOutput("gzip -t ".$fileToTest) eq "0"){
    return 1;
  }else{
    return 0;
  }
}


my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $adptsim = $pathdir."/src/adptSim";
my $deamsim = $pathdir."/src/deamSim";
my $fragsim = $pathdir."/src/fragSim";
my $artprog = $pathdir."/art_src_MountRainier/art_illumina";





fileExists($adptsim);
fileExists($deamsim);
fileExists($fragsim);
fileExists($artprog);

my $fraglength=35;
my $numberOfFragments=1000;
my $comp = "0,0,1";
my $minsize=0;
my $maxsize=1000;
my $adapterF="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
my $adapterR="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
my $readlength=75;
my $distmis=1;

my $qs =0;
my $qs2=0;


sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is a wrapper to run the different programs to create
 a set of Illumina reads for ancient DNA from a set of fasta files
 representing the endogenous, the contaminant from the same species
 and the bacterial contamination.

\n\n usage:\t".$0." <options> [directory with fasta directories] \n\n".
" This directory should contain 3 directories:\n".
" \tbact/ The bacterial contamination\n".
" \tcont/ The contamination from the same species\n".
" \tendo/ The endogenous material\n".
" \n".

" Options:\n".

  "\t--comp [B,C,E]\t\t\t\tComposition of the final set in fraction \n".
  "\t\t\t\t\t\tthe 3 numbers represent the bacterial, contaminant and endogenous\n".
  "\t\t\t\t\t\tex: --comp 0.6,0.02,0.38 will result\n".
  "\t\t\t\t\t\tin 60% bacterial contamination while the rest will be from the same\n".
  "\t\t\t\t\t\tspecies 5% will be contamination and 95% will be endogenous\n".
  "\t\t\t\t\t\tDefault: --comp ".$comp."\n".
 "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".
  "\t-o\t\t\t\t\tOutput prefix (default: [input dir]/simadna)\n".
  " Either specify:\n".
  "\t\t-n\t[number]\t\tGenerate [number] fragments (default: ".$numberOfFragments.")\n".
  "\t\t-c\t[coverage]\t\tEndogenous coverage\n".
  "\n".
  " Fragment selection\n".
  " ===================\n".

  " \t\t--misince\t[file]\t\tBase misincorporation for the endogenous fragments (default none)\n".
  " \t\t--misincc\t[file]\t\tBase misincorporation for the contaminant fragments (default none)\n".
  " \t\t--misincb\t[file]\t\tBase misincorporation for the bacterial fragments (default none)\n".
  " \t\t--distmis\t[dist]\t\tDistance to consider for base misincorporation (default ".$distmis.")\n".

  " \t\t\t\t\t\tthis file is obtained from mapDamage\n".
  "\n".
  " \tFragment size distribution: specify either one of the 3 possible options:
		-l\t[length]\t\tGenerate fragments of fixed length  (default: ".$fraglength.")
		-s\t[file]\t\t\tOpen file with size distribution (one fragment length per line)
		-f\t[file]\t\t\tOpen file with size frequency in the following format:
				     \t\tlength[TAB]freq	ex:
				     \t\t40	0.0525
				     \t\t41	0.0491
				     \t\t...

		Length options:
		\t--loc\t[file]\t\tLocation for lognormal distribution (default none)
		\t--scale\t[file]\t\tScale for lognormal distribution    (default none)\n\n".
  " \tFragment size limit:
		--minsize\t[length]\tMinimum fragments length (default: ".$minsize.")
                --maxsize\t[length]\tMaximum fragments length (default: ".$maxsize.")\n".
  "\n".
  " \tFragment methylation:
		--methyl\t\t\tAllow for lowercase C and G to denote
                \t\t\t\tmethylated cytosines on the + and - strand respectively
                \t\t\t\t(default: not used)\n".
  "\n".
  " Deamination\n".
  " ===================\n".
  " To add deamination to the bacterial and endogenous material,you can specify\n".
  " either one of these options:\n".

  "\t-mapdamage\t[mis.txt] [protocol]\tRead the miscorporation file [mis.txt]
                           	       \t\tproduced by mapDamage
		                       \t\t[protocol] can be either \"single\" or \"double\" (without quotes)
		                       \t\tSingle strand will have C->T damage on both ends
		                       \t\tDouble strand will have and C->T at the 5' end and G->A damage at the 3' end".
  "\n".
  "\t-matfile\t[matrix file prefix]\tRead the matrix file of substitutions
		                                Provide the prefix only, both files must end with
		                               	5.dat and 3.dat\n".


  "\t-damage\t\t[v,l,d,s]	  \tFor the Briggs et al. 2007 model
		                  \t\tThe parameters must be comma-separated e.g.: -damage 0.03,0.4,0.01,0.3
		                  \t\t\tv: nick frequency
		                  \t\t\tl: length of overhanging ends (geometric parameter)
		                  \t\t\td: prob. of deamination of Cs in double-stranded parts
		                  \t\t\ts: prob. of deamination of Cs in single-stranded parts\n".


  "\n".
  " Alternatively, you can specify these options independently for the endogenous (e), bacterial (b)\n".
  " and present-day human contaminant (c) using the following options:\n".

  "\t-mapdamagee\t[mis.txt] [protocol]\tEndogenous mapDamage misincorporation file\n".
  "\t-matfilee\t[matrix file prefix]\tEndogenous matrix file of substitutions\n".
  "\t-damagee\t[v,l,d,s]	  \tEndogenous Briggs parameters\n".

  "\t-mapdamageb\t[mis.txt] [protocol]\tBacterial mapDamage misincorporation file\n".
  "\t-matfileb\t[matrix file prefix]\tBacterial matrix file of substitutions\n".
  "\t-damageb\t[v,l,d,s]	  \tBacterial Briggs parameters\n".

  "\t-mapdamagec\t[mis.txt] [protocol]\tHuman contaminant mapDamage misincorporation file\n".
  "\t-matfilec\t[matrix file prefix]\tHuman contaminant matrix file of substitutions\n".
  "\t-damagecd\t[v,l,d,s]	  \tHuman contaminant Briggs parameters\n".

  "\n please note that if you do specify deamination for one source but not for another, no deamination will be added\n".

  "\n".
  " If using --methyl, you can also specify different matrix file for methylated\n".
  "\t-matfilenonmeth\t[matrix file prefix]\tRead the matrix file of substitutions for non-methylated Cs
		                            \tProvide the prefix only, both files must end with
		                            \t5.dat and 3.dat
   \t-matfilemeth\t[matrix file prefix]\tRead the matrix file of substitutions for methylated Cs
		                         \tProvide the prefix only, both files must end with
		                         \t5.dat and 3.data".

  "\n".
  " Adapter and sequencing\n".
  " ===================\n".
  "	-fa	[seq]			\tAdapter that is observed after the forward read (Default: ".substr($adapterF,0,10)."...)
	-sa	[seq]			\tAdapter that is observed after the reverse read (Default: ".substr($adapterR,0,10)."...)
	-rl	[length]		\tDesired read length  (Default: ".$readlength.")
	-se                             \tuse single-end sequencing (Default: paired-end)

                                        \tThe following options change the sequencing error rate, please note that positive factor
                                        \twill decrease the rate of such errors and a negative one will increase it.
        -qs     [factor]                \tIncrease error rate for forward reads by a factor of 1/(10^([factor]/10)) (Default: ".$qs.")
        -qs2    [factor]                \tIncrease error rate for reverse reads by a factor of 1/(10^([factor]/10)) (Default: ".$qs2.")


	-ss     [system]                \tIllumina platfrom to use, the parentheses indicate the max. read length
	                                \tuse the shorthand in the left column:
                                                                        (single-end, paired-end)
\t\t\t\t\t\t   GA2  - GenomeAnalyzer II (  50bp,  75bp)
\t\t\t\t\t\t   HS20 - HiSeq 2000        ( 100bp, 100bp)
\t\t\t\t\t\t   HS25 - HiSeq 2500        ( 125bp, 150bp) (Default)
\t\t\t\t\t\t   HSXt - HiSeqX TruSeq     ( 150bp, 150bp)
\t\t\t\t\t\t   MSv1 - MiSeq v1          ( 250bp, 250bp)
\t\t\t\t\t\t   MSv3 - MiSeq v3          ( 250bp, 250bp)".

	  "\n\n".
  "\n\n".
#"  Output Options:\n".
#"  Input Options:\n".
"\n\n";
  exit;
}

my $uniq =0;

my $starttime = time;

my $coverage=undef;

my $help;
my $outputprefix;
my $compB;
my $compC;
my $compE;

my $filefragsize;
my $filefragfreqsize;
my $loc;
my $scale;
my $misince;
my $misincb;
my $misincc;


my @mapdamage;
my $matfile;
my $matfilenonmeth;
my $matfilemeth;

my $briggs;

my @mapdamagee;
my $matfilee;
my $briggse;

my @mapdamageb;
my $matfileb;
my $briggsb;

my @mapdamagec;
my $matfilec;
my $briggsc;



my $se=0;
my $ss;

my $fa;
my $sa;
my $rl;

my $methyl;
my $dirWithChr = $ARGV[$#ARGV];

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'mock' => \$mock, 'methyl' => \$methyl, 'uniq' => \$uniq, 'se' => \$se, 'ss=s' => \$ss, 'distmis=i' => \$distmis, 'misince=s' => \$misince,'misincb=s' => \$misincb,'misincc=s' => \$misincc, 'comp=s' => \$comp,'mapdamage=s{2}' => \@mapdamage, 'mapdamagee=s{2}' => \@mapdamagee, 'mapdamageb=s{2}' => \@mapdamageb, 'mapdamagec=s{2}' => \@mapdamagec,'matfile=s' => \$matfile, 'damage=s' => \$briggs,'matfilee=s' => \$matfilee,'matfilenonmeth=s' => \$matfilenonmeth, ,'matfilemeth=s' => \$matfilemeth, 'damagee=s' => \$briggse,'matfileb=s' => \$matfileb, 'damageb=s' => \$briggsb,'matfilec=s' => \$matfilec, 'damagec=s' => \$briggsc,'o=s' => \$outputprefix, 'n=i' => \$numberOfFragments,'l=i' => \$fraglength, 's=s' => \$filefragsize, 'f=s' => \$filefragfreqsize, 'loc=s' => \$loc, 'fa=s' => \$fa, 'sa=s' => \$sa, 'rl=s' => \$rl, 'scale=s' => \$scale, 'c=f' => \$coverage, 'minsize=i' => \$minsize,'maxsize=i' => \$maxsize,'qs=i' => \$qs,'qs2=i' => \$qs2)
          or defined $help );

if( !(defined $ss) ){
  $ss = "HS25";
}

if( defined $fa ){
  $adapterF=uc($fa);
}

if( defined $sa ){
  $adapterR=uc($sa);
}

if( defined $rl ){
  $readlength=$rl;
}

if(index($adapterF, "N") != -1) {
  my $tempadapStr="";
  my @DNAalphabet = ("A","C","G","T");

  for(my $i=0;$i<length($adapterF);$i++){
    if(substr($adapterF,$i,1) eq "N"){
      $tempadapStr=$tempadapStr.$DNAalphabet[int(rand(4))];
    }else{
      $tempadapStr=$tempadapStr.substr($adapterF,$i,1);
    }
  }
  warn "WARNING: Unresolved bases found in forward adapter\n"."Replaced:\n".$adapterF."\nby random bases:\n".$tempadapStr."\n";
  $adapterF=$tempadapStr;
}

if(index($adapterR, "N") != -1) {
  my $tempadapStr="";
  my @DNAalphabet = ("A","C","G","T");

  for(my $i=0;$i<length($adapterR);$i++){
    if(substr($adapterR,$i,1) eq "N"){
      $tempadapStr=$tempadapStr.$DNAalphabet[int(rand(4))];
    }else{
      $tempadapStr=$tempadapStr.substr($adapterR,$i,1);
    }
  }
  warn "WARNING Unresolved bases found in reverse adapter\n"."Replaced:\n".$adapterR."\nby random bases:\n".$tempadapStr."\n";
  $adapterR=$tempadapStr;
}




if ($ss eq "GA2") {		#- GenomeAnalyzer II (50bp, 75bp)
  if ($se) {
    if ($readlength>50) {
      die "Read length ".$readlength." is greater than the one allowed by the platform\n";
    }
  } else {
    if ($readlength>75) {
      die "Read length ".$readlength." is greater than the one allowed by the platform\n";
    }
  }
} else {
  if ($ss eq "HS20") {		#- HiSeq 2000 (100bp,100bp)
    #if ($se) {
      if ($readlength>100) {
	die "Read length ".$readlength." is greater than the one allowed by the platform\n";
      }
    #} else {
    #  die "The platform does not provide paired-end sequencing\n";
    #}
  } else {
    if ($ss eq "HS25") {	#- HiSeq 2500 (125bp, 150bp) (Default)
      if ($se) {
	if ($readlength>125) {
	  die "Read length ".$readlength." is greater than the one allowed by the platform\n";
	}
      } else {
	if ($readlength>150) {
	  die "Read length ".$readlength." is greater than the one allowed by the platform\n";
	}
      }
    } else {
      if ($ss eq "HSXt") {	#- HiSeqX TruSeq (150bp,150bp)
	#$if ($se) {
	if ($readlength>150) {
	  die "Read length ".$readlength." is greater than the one allowed by the platform\n";
	}
	#} else {
	#  die "The platform does not provide paired-end sequencing\n";
	#}
      } else{
	if ($ss eq "MSv1") {	#- MiSeq v1 (250bp,250bp)
	  #if ($se) {
	  if ($readlength>250) {
	    die "Read length ".$readlength." is greater than the one allowed by the platform\n";
	  }
	#} else {
	#    die "The platform does not provide paired-end sequencing\n";
	  #  }
	} else {
	  if ($ss eq "MSv3") {	#- MiSeq v3 (250bp,250bp)
	    #if ($se) {
	    if ($readlength>250) {
	      die "Read length ".$readlength." is greater than the one allowed by the platform\n";
	    }
	    #} else {
	    #  die "The platform does not provide paired-end sequencing\n";
	    #}
	  } else {
	    die "Invalid sequencing platform ".$ss."\n";
	  }
	}#not msv1
      }#not HSXt
    }#not HS25
  }#not HS20
}#not GA2

if( (defined $misince) ){
  if( !(-e $misince)){
    die "Endogenous misincorporation file does not exists\n";
  }
}

if( (defined $misincc) ){
  if( !(-e $misincc)){
    die "Contaminant misincorporation file does not exists\n";
  }
}

if( (defined $misincb) ){
  if( !(-e $misincb)){
    die "Bacterial misincorporation file does not exists\n";
  }
}

if( ($qs < -93) || ($qs > 93) ){
  die "Please enter an integer between -93 and 93 for qs\n";
}

if( ($qs2 < -93) || ($qs2 > 93) ){
  die "Please enter an integer between -93 and 93 for qs2\n";
}

if ($se) {
  if( ($qs2 != 0 ) ){
    die "Do not enter a value for -qs2 when using single-end\n";
  }
}

checkDeamParam($briggs, $matfile );
checkDeamParam($briggse,$matfilee);
checkDeamParam($briggsb,$matfileb);
checkDeamParam($briggsc,$matfilec);



if ( (defined $matfile          ) ||
     (defined $matfilenonmeth   ) ||
     (defined $matfilemeth      ) ||
     (defined $briggs           ) ||
     (@mapdamage ) ) {

  if ( (defined $matfilee) ||
       (defined $briggse )  ) {
    die "Specify either patterns for endogenous and bacterial using either -matfile or -damage but do not specify other parameters\n";
  }

  if ( (defined $matfileb) ||
       (defined $briggsb )  ) {
    die "Specify either patterns for endogenous and bacterial using either -matfile or -damage but do not specify other parameters\n";
  }

  if ( (defined $matfilec) ||
       (defined $briggsc )  ) {
    die "Specify either patterns for endogenous and bacterial using either -matfile or -damage but do not specify other parameters\n";
  }


  if(defined $matfile){
    $matfilee = $matfile;
    $matfileb = $matfile;
  }

  if(defined $briggs){
    $briggse = $briggs;
    $briggsb = $briggs;
  }

  if(defined $matfilenonmeth){
    if( not(defined $methyl) ){
      die "Must define --methyl for the fragment selection\n";
    }
  }

  if(defined $matfilemeth){
    if( not(defined $methyl) ){
      die "Must define --methyl for the fragment selection\n";
    }
  }

  if( (not(defined $matfilenonmeth) &&
       (   defined $matfilemeth) )
      ||
      (   (defined $matfilenonmeth) &&
       not(defined $matfilemeth) )){
    die "Must defined matrix file for both methylated and non-methylated";
  }

  if( (defined $matfilenonmeth) ||
      (defined $matfilemeth   )  ){

    if( (defined $matfilee) ||
	(defined $matfileb) ||
	(defined $matfilec)  ) {
      die "Cannot specify as of yet, different deamination profiles for the endogenous, bacterial, contaminant and methylation deamination\n";
    }

  }


  if(@mapdamage){
    @mapdamagee = @mapdamage;
    @mapdamageb = @mapdamage;
  }



}

if( defined $methyl ){
  if( not(defined $matfilenonmeth) ){
    die "Must define --matfilenonmeth if methylation is used\n";
  }
}

if( defined $methyl ){
  if( not(defined $matfilemeth) ){
    die "Must define --matfilemeth if methylation is used\n";
  }
}


















if( (defined $loc) &&
    !(defined $scale) ){
  die "Must specify both --loc and --scale, not just one";
}

if( !(defined $loc) &&
    (defined $scale) ){
  die "Must specify both --loc and --scale, not just one";
}

if( (defined $filefragsize) ){
  if($fraglength!=35){
    die "Cannot specify both -l and -s";
  }
}

if( (defined $filefragfreqsize) ){
  if($fraglength!=35){
    die "Cannot specify both -l and -f";
  }
}

if( (defined $filefragsize)     &&
    (defined $filefragfreqsize) ){
  die "Cannot specify both -s and -f";
}


my @arraycomp = split(",",$comp);
if($#arraycomp != 2){
  die "the --comp option must have 3 comma-delimited fields, found ".($#arraycomp+1);
}

if( (defined $coverage) &&
    $numberOfFragments!=1000){
  die "Must use the -c or -n but not both";
}



$compB=$arraycomp[0];
$compC=$arraycomp[1];
$compE=$arraycomp[2];
if(!looks_like_number($compB)){
  die "The --comp option must have 3 comma-delimited numbers";
}
if(!looks_like_number($compC)){
  die "The --comp option must have 3 comma-delimited numbers";
}
if(!looks_like_number($compE)){
  die "The --comp option must have 3 comma-delimited numbers";
}

if( ($compB+$compC+$compE) != 1){
  die "The --comp option must have 3 comma-delimited numbers that sum up to 1, current sum ".($compB+$compC+$compE);
}


if( (defined $coverage) && ($compE==0) ){
  die "Cannot use the endogenous coverage -c and not specify some level of endogenous content.";
}


if(substr($dirWithChr,length($dirWithChr)-1,1) ne "/"){
  $dirWithChr = $dirWithChr."/";
}

if( !(defined $outputprefix) ){
  $outputprefix = $dirWithChr."simadna";
}else{
}


my @arrayofdirs=listdir($dirWithChr);

@arrayofdirs=sort(@arrayofdirs);
if($#arrayofdirs != 2){
  die "The input directory must contain 3 directories:
 bact/
 endo/
 cont/
";
}

if($arrayofdirs[0] ne $dirWithChr."bact/" ||
   $arrayofdirs[1] ne $dirWithChr."cont/" ||
   $arrayofdirs[2] ne $dirWithChr."endo/" ){
  die "The input directory must contain 3 directories named:
 bact/
 endo/
 cont/
";
}

my @arrayofFilesbact = listdirFa(   $arrayofdirs[0] );
my @arrayofFilescont = listdirFa(   $arrayofdirs[1] );
my @arrayofFilesendo = listdirFa(   $arrayofdirs[2] );
my @arrayofCdirsendo = listdirCdir( $arrayofdirs[2] );

my $sumB=0;
my $sumC=0;
my $sumE=0;

my @arrayofFilesbactSL;
my @arrayofFilesbactL;
my @arrayofFilesbactLfrac;
my @arrayofFilesbactToExtract;

my @arrayofFilesbactLFromList;
my @arrayofFilesbactLFromListW;
my @arrayofFilesbactLFromListB;
my @arrayofFilesbactLFromListP;

if($compB>0){

  if($#arrayofFilesbact==-1){
    die "If you want bacterial contamination, please have at least one file in the bact/ directory\n";
  }

  if( !( -e $arrayofdirs[0]."/list" ) ){
    die "List file ".$arrayofdirs[0]."/list does not exist, this file must contain the list of files and their weight ex:
file1.fa\t0.3
file2.fa\t0.2
file3.fa\t0.15
file4.fa\t0.12
file5.fa\t0.1
file6.fa\t0.7
file7.fa\t0.5
\n";
  }else{
    my $sumWeight=0;
    open(FILE,$arrayofdirs[0]."/list");
    while(my $line = <FILE>){
      if($line =~ /^(\S+)\s+(\S+)$/){

	if(!looks_like_number($2)){
	  die "Cannot parse line from bacterial list: ".$arrayofdirs[0]."/list must be:\nfile\tweight(between 0 and 1)\nfound line:".$line;
	}
	push(@arrayofFilesbactLFromList, $1);
	push(@arrayofFilesbactLFromListB, 0);
	push(@arrayofFilesbactLFromListW,$2);
	$sumWeight+=$2;
      }else{
	die "Cannot parse line from bacterial list: ".$arrayofdirs[0]."/list line:".$line;
      }
    }
    close(FILE);

    if($sumWeight<0.99 || $sumWeight>1.01 ){
      die "Problem from bacterial list: ".$arrayofdirs[0]."/list sum is not 1, found: ".$sumWeight;
    }

  }

}

foreach my $fafile (@arrayofFilesbact){
  print STDERR "Found bacterial contaminant file ".$fafile."\n";

  if(isZipped($fafile)){
    die "The following file ".$fafile." is zipped, please unzip it\n";
  }

  if(!(-f $fafile.".fai")){
    my $cmd = "samtools faidx $fafile";
    runcmdforce($cmd);
  }else{
    print STDERR "Found bacterial contaminant indx ".$fafile.".fai\n";
  }
  my $sumForFile=0;
  open(FILE,$fafile.".fai");
  while(my $line = <FILE>){
    if($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
      $sumB+=$2;
      $sumForFile+=$2;
    }else{
      die "Cannot parse line from fasta index ".$fafile.".fai";
    }
  }
  close(FILE);

  for(my $i=0;$i<=$#arrayofFilesbactLFromList;$i++){
    print $arrayofdirs[0].$arrayofFilesbactLFromList[$i]."\n";
    print $fafile."\n";
    if($arrayofdirs[0].$arrayofFilesbactLFromList[$i] eq $fafile){
      if($arrayofFilesbactLFromListB[$i] == 1){
	die "Fasta  ".$fafile." was found twice in the list";
      }
      $arrayofFilesbactLFromListB[$i] = 1;
    }
  }


  push(@arrayofFilesbactL, $sumForFile);
  push(@arrayofFilesbactSL,$sumB);
  push(@arrayofFilesbactToExtract,0);
}


for(my $i=0;$i<=$#arrayofFilesbactLFromList;$i++){
  if($arrayofFilesbactLFromListB[$i] == 0){
    die "Fasta  ".$arrayofFilesbactLFromList[$i]." was not found in the directory but it was found in the list\n";
  }
}


my $sumOfListW=0;
push(@arrayofFilesbactLFromListP,$sumOfListW);

for(my $i=0;$i<=$#arrayofFilesbactLFromListW;$i++){
  #print $i."\t".$arrayofFilesbactLFromListW[$i]."\n";
  $sumOfListW += $arrayofFilesbactLFromListW[$i];
  push(@arrayofFilesbactLFromListP,$sumOfListW);
}




if($compC>0){
  if($#arrayofFilescont==-1){
    die "If you want contamination, please have at least one file in the cont/ directory\n";
  }
}

my @arrayofFilescontSL;
my @arrayofFilescontL;
my @arrayofFilescontLfrac;
my @arrayofFilescontToExtract;

foreach my $fafile (@arrayofFilescont){
  print STDERR "Found present-day human contamination file  ".$fafile."\n";


  if(isZipped($fafile)){
    die "The following file ".$fafile." is zipped, please unzip it\n";
  }


  if(!(-f $fafile.".fai")){
    my $cmd = "samtools faidx $fafile";
    runcmdforce($cmd);
  }else{
    print STDERR "Found present-day human contamination index ".$fafile.".fai\n";
  }
  my $sumForFile=0;
  open(FILE,$fafile.".fai");
  while(my $line = <FILE>){
    if($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
      $sumC+=$2;
      $sumForFile+=$2;
    }else{
      die "Cannot parse line from fasta index ".$fafile.".fai";
    }
  }
  close(FILE);
  push(@arrayofFilescontL, $sumForFile);
  push(@arrayofFilescontSL,$sumC);
  push(@arrayofFilescontToExtract,0);
}

foreach my $s (@arrayofFilescontL){
  push( @arrayofFilescontLfrac, ($s/$sumC) );
}

my $diploidMode=0;

my $multiCellMode   =0;
my $multiCellModeMAX=0;

my @arrayOfEndoCellFasta;

if ($compE>0) {			#if we have endogenous material

  #my @arrayofFilesendo = listdirFa(   $arrayofdirs[2] );
  #my @arrayofCdirsendo = listdirCdir( $arrayofdirs[2] );
  my $cell0Dir;

  #
  #  MULTICELL MODE
  #
  if ( ($#arrayofCdirsendo) >= 0 ) {
    $multiCellMode = 1;
    if ( $#arrayofFilesendo != -1) {
      die "The endogenous directory ".$arrayofdirs[2]." must directories called CX where X is 0...N without any other fasta files, found: ".($#arrayofFilesendo+1)." fasta files\n";
    }


    my $foundC0=0;
    my $foundCMAX=0;

    foreach my $dirCendo (@arrayofCdirsendo){

      if($dirCendo =~ /C(\d+)$/        ){
	if($1 > $foundCMAX){
	  $foundCMAX = $1;
	}
      }else{
	die "Directories ".$arrayofdirs[2]." does not correspond to the pattern CX where X is 0...N, found: $dirCendo\n";
      }

      if($dirCendo =~ /C0$/){
	$foundC0=1;
      }

    }

    if(!$foundC0){
      die "Directories ".$arrayofdirs[2]." does not correspond to the pattern CX where X is 0...N, must start with C0\n";
    }

    $multiCellModeMAX=$foundCMAX;

    for(my $cendoI=0;$cendoI<=$foundCMAX;$cendoI++){
      if( grep { $_ eq $arrayofdirs[2]."C".$cendoI."" } @arrayofCdirsendo ){
	#fine
      }else{
	die "Missing directory ".$arrayofdirs[2]."C".$cendoI."\n";
      }

      my @arrayofFilesendoCell = listdirFa(  $arrayofdirs[2]."C".$cendoI."/");
      $arrayOfEndoCellFasta[ $cendoI ] = \@arrayofFilesendoCell;

      if ($cendoI == 0 ) {
	@arrayofFilesendo = @arrayofFilesendoCell;
	$cell0Dir = $arrayofdirs[2]."C".$cendoI."/";
	if ( ($#arrayofFilesendoCell+1) == 1) {
	  $diploidMode=0;
	} else {
	  if ( ($#arrayofFilesendoCell+1) == 2) {
	    $diploidMode=1;
	  } else {
	    die "The endogenous directory must have 1 (haploid) or 2 (diploid) files, directory in question: ".($arrayofdirs[2]."C".$cendoI).", found ".($#arrayofFilesendoCell+1)."\n";
	  }
	}
      }else{
	if($diploidMode == 0 ){
	  if ( ($#arrayofFilesendoCell+1) != 1) {
	    die "The endogenous directory must have all have 1 (haploid) or 2 (diploid) files, directory in question: ".($arrayofdirs[2]."C".$cendoI).", found: ".($#arrayofFilesendoCell+1)."\n";
	  }
	}

	if($diploidMode == 1 ){
	  if ( ($#arrayofFilesendoCell+1) != 2) {
	    die "The endogenous directory must have all have 1 (haploid) or 2 (diploid) files, directory in question: ".($arrayofdirs[2]."C".$cendoI).", found: ".($#arrayofFilesendoCell+1)."\n";
	  }
	}


      }
    }

    print STDERR  "Found ".($#arrayofCdirsendo+1)." directories: C[0..C".$foundCMAX."]\n";


  } else {
    #
    #  SINGLE CELL MODE
    #

    if ( ($#arrayofFilesendo+1) == 1) {
      $diploidMode=0;
    } else {
      if ( ($#arrayofFilesendo+1) == 2) {
	$diploidMode=1;
      } else {
	die "The endogenous directory must have 1 (haploid) or 2 (diploid) files, found: ".($#arrayofFilesendo+1)."";
      }
    }

  }

  foreach my $fafile (@arrayofFilesendo) {
    print STDERR "Found endogenous file  ".$fafile."\n";
    if (!(-f $fafile.".fai")) {
      my $cmd = "samtools faidx $fafile";
      runcmdforce($cmd);
    } else {
      print STDERR "Found endogenous index ".$fafile.".fai\n";
    }

    open(FILE,$fafile.".fai");
    while (my $line = <FILE>) {
      if ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	$sumE+=$2;
      } else {
	die "Cannot parse line from fasta index ".$fafile.".fai";
      }
    }
    close(FILE);
  }


  #comparing endogenous fai
  my @arrayofFilesendofai;
  if($multiCellMode){
    @arrayofFilesendofai = listdirFai( $cell0Dir );
  }else{
    @arrayofFilesendofai = listdirFai( $arrayofdirs[2] );
  }

  if ($diploidMode) {
    if ( ($#arrayofFilesendofai+1) != 2) {
      die "The endogenous directory must have 2 fai files, found: ".($#arrayofFilesendofai+1)." files";
    }

    my @arrayFaiChr1;
    my @arrayFaiChr2;

    open(FILE,$arrayofFilesendofai[0]);
    while (my $line = <FILE>) {
      #              1       2       3       4       5
      if ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	push(@arrayFaiChr1,$1);
      } else {
	die "Cannot parse line from fasta index ".$arrayofFilesendofai[0].".fai";
      }
    }
    close(FILE);

    open(FILE,$arrayofFilesendofai[1]);
    while (my $line = <FILE>) {
      #              1       2       3       4       5
      if ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	push(@arrayFaiChr2,$1);
      } else {
	die "Cannot parse line from fasta index ".$arrayofFilesendofai[1].".fai";
      }
    }
    close(FILE);

    if ($#arrayFaiChr1 != $#arrayFaiChr2 ) {
      die "The endogenous directory must have 2 fai files, with the same chromosomes, one has ".$#arrayFaiChr1." while the other has ".$#arrayFaiChr2."\n";
    }

    for (my $i=0;$i<=$#arrayFaiChr1;$i++) {
      if ( $arrayFaiChr1[$i] ne $arrayFaiChr2[$i] ) {
	die "The endogenous directory must have 2 fai files, with the same chromosomes, one has ".$arrayFaiChr1[$i]." while the other has ".$arrayFaiChr2[$i]."\n";

      }
    }

  }


}

print STDERR "\nFound ".$sumB." bp bacterial, ".$sumC." bp contaminant, ".$sumE." bp endogenous\n";

my $averageSize;

if(defined $filefragfreqsize){
  my $num=0;
  my $sum=0;
  open(FILE,"gunzip -f -c ".$filefragfreqsize." | ");
  while(my $line = <FILE>){
    chomp($line);
    my @arr = split("\t",$line);
    if($#arr != 1){
      die "Line $line in $filefragfreqsize should have 2 tab-delimited fields\n";
    }
    $sum+=$arr[0]*$arr[1];
    $num+=$arr[1];
  }
  close(FILE);
  $averageSize=($sum/$num);
} else {
  if (defined $filefragsize) {
    my $num=0;
    my $sum=0;
    open(FILE,"gunzip -f -c ".$filefragsize." | ");
    while (my $line = <FILE>) {
      chomp($line);
      my @arr = split("\t",$line);
      if ($#arr != 0) {
	die "Line $line in file $filefragsize should have 1 tab-delimited fields\n";
      }

      $sum+=$line;
      $num++;
    }
    close(FILE);
    $averageSize=($sum/$num);

  } else {
    if (defined $loc) {
      $averageSize = exp( $loc+($scale**2)/2);
    } else {
      $averageSize = $fraglength;
    }
  }
}

print STDERR "Computed average size: ".$averageSize."bp\n";

my $numberOfFragmentsC;
my $numberOfFragmentsB;
my $numberOfFragmentsE;

if(defined $coverage){
  my $sumEeffective=$sumE;
  if($diploidMode){
    $sumEeffective=$sumE/2;#if diploid, the sum is overestimated by 2.
  }
  $numberOfFragmentsE = int($coverage*($sumEeffective/$averageSize));
  $numberOfFragments  = int($numberOfFragmentsE/$compE);
  $numberOfFragmentsC = int($numberOfFragments * $compC/($compE+$compC+$compB));
  $numberOfFragmentsB = int($numberOfFragments * $compB/($compE+$compC+$compB));
}else{
  $numberOfFragmentsC = int($numberOfFragments*$compC);
  $numberOfFragmentsB = int($numberOfFragments*$compB);
  $numberOfFragmentsE = int($numberOfFragments*$compE);
}

$numberOfFragments =   $numberOfFragmentsE+$numberOfFragmentsC+$numberOfFragmentsB;


my $numberOfFragmentsE1=0;
my $numberOfFragmentsE2=0;

if($diploidMode){
  for(my $i=0;$i<$numberOfFragmentsE;$i++){
    if(rand()<0.5){
      $numberOfFragmentsE1++;
    }else{
      $numberOfFragmentsE2++;
    }
  }
}else{
  $numberOfFragmentsE1 = $numberOfFragmentsE;
  $numberOfFragmentsE2 = 0;
}








if($numberOfFragmentsB>0){
  print STDERR "\nSelecting bacterial genomes according to composition vector:\n";
  my $maxsizeFilename=0;

  for(my $i=0;$i<=$#arrayofFilesbactLFromList;$i++){
    if(length($arrayofFilesbactLFromList[$i])>$maxsizeFilename){
      $maxsizeFilename = length($arrayofFilesbactLFromList[$i]);
    }
  }
  for(my $i=0;$i<=$#arrayofFilesbactLFromList;$i++){
    print STDERR "P[".(' 'x($maxsizeFilename-length($arrayofFilesbactLFromList[$i]))).$arrayofFilesbactLFromList[$i]."] = ".$arrayofFilesbactLFromListW[$i]."\n";
  }

}

if (0) {
  for (my $i=0;$i<$numberOfFragmentsB;$i++) {
    my $randB=int(rand($sumB));

    for (my $j=0;$j<=$#arrayofFilesbactSL;$j++) {
      if ($randB<=$arrayofFilesbactSL[$j]) {
	$arrayofFilesbactToExtract[$j]++;
	last;
      }
    }
  }
} else {

  if ($numberOfFragmentsB>0) {
    my $digitfb= log($numberOfFragmentsB)/log(10) +1;
    for (my $i=0;$i<$numberOfFragmentsB;$i++) {
      if ($i!=0&&
	  ($i%100000)==0) {
	print  STDERR "selected ".sprintf("%".$digitfb."s",$i)." bacterial fragments\n";
      }

      #print  $i."\n";
      my $randP=rand();
      my $indexFound=-1;
      #print  $randP."\n";
      for (my $i=1;$i<=$#arrayofFilesbactLFromListP;$i++) {
	#print "prob ".$randP." ".$arrayofFilesbactLFromListP[$i-1]." ".$arrayofFilesbactLFromListP[$i];
	if ( ($arrayofFilesbactLFromListP[$i-1]<=$randP) &&
	     ($arrayofFilesbactLFromListP[$i]  >=$randP) ) {
	  $indexFound=$i-1;
	  #die;
	  last;
	}
      }

      if ($indexFound != -1) {
	$arrayofFilesbactToExtract[$indexFound]++;
      } else {
	$i--;			#restart iteration
	#die;
      }
    }
  }

}


print STDERR "We will generate: \n\n";
print STDERR "".$numberOfFragments."\ttotal fragments\n";
print STDERR "--------------------------------------------\n";
print STDERR "".$numberOfFragmentsE."\t(".sprintf("% .2f",100*$numberOfFragmentsE/$numberOfFragments)."%) "."\tendogenous fragments\n";

my @arrayFragForEachCell;
my @arrayFragForEachCellE1;
my @arrayFragForEachCellE2;

if($multiCellMode){

  for(my $i=0;$i<=$multiCellModeMAX;$i++){
    $arrayFragForEachCell[$i] = 0;
    $arrayFragForEachCellE1[$i] = 0;
    $arrayFragForEachCellE2[$i] = 0;
  }

  #force the first file to have a least a fragment
  if($numberOfFragmentsE>0){
    $arrayFragForEachCell[0]++;
  }

  for(my $i=0;$i<$numberOfFragmentsE;$i++){
    my $randB=int(rand($multiCellModeMAX+1));
    $arrayFragForEachCell[$randB]++;
  }

  my $digE = int(log($numberOfFragmentsE)/log(10))+2;

  if($numberOfFragmentsE>0){
    $arrayFragForEachCellE1[0] ++;
  }

  for(my $i=0;$i<=$multiCellModeMAX;$i++){
    if($diploidMode){
      for(my $j=0;$j<$arrayFragForEachCell[$i];$j++){
	if(rand()<0.5){
	  $arrayFragForEachCellE1[$i]++;
	}else{
	  $arrayFragForEachCellE2[$i]++;
	}
      }
    }else{
      $arrayFragForEachCellE1[$i] = $arrayFragForEachCell[$i];
      $arrayFragForEachCellE2[$i] = 0;
    }


    print STDERR "Will extract ".sprintf("%".$digE."d",$arrayFragForEachCell[$i])." from cell: ".$arrayofdirs[2]."C".$i."/\t(".sprintf("% .2f",100*$arrayFragForEachCell[$i]/$numberOfFragmentsE)."%) "."\tendogenous fragments\n";
  }
  #die;
} else {

  if ($diploidMode) {
    print STDERR "".$numberOfFragmentsE1."\t(".sprintf("% .2f",100*$numberOfFragmentsE1/$numberOfFragments)."%) "."\tendogenous fragments from first  chr file: ".$arrayofFilesendo[0]."\n";
    print STDERR "".$numberOfFragmentsE2."\t(".sprintf("% .2f",100*$numberOfFragmentsE2/$numberOfFragments)."%) "."\tendogenous fragments from second chr file: ".$arrayofFilesendo[1]."\n";
  } else {
    #nothing to print
  }

}


print STDERR "--------------------------------------------\n";
print STDERR   "".$numberOfFragmentsC."\t(".sprintf("% .2f",100*$numberOfFragmentsC/$numberOfFragments)."%) "."\tcontaminant fragments\n";


for(my $i=0;$i<$numberOfFragmentsC;$i++){
  my $randC=int(rand($sumC));

  for(my $j=0;$j<=$#arrayofFilescontSL;$j++){
    if($randC<=$arrayofFilescontSL[$j]){
      $arrayofFilescontToExtract[$j]++;
      last;
    }
  }

}

for(my $i=0;$i<=$#arrayofFilescont;$i++){
  print STDERR "".$arrayofFilescontToExtract[$i]."\t(".sprintf("% .2f",100*$arrayofFilescontToExtract[$i]/$numberOfFragments)."%) "."\tcontaminant fragments from file: ".$arrayofFilescont[$i]."\n";
}

print STDERR "--------------------------------------------\n";
print STDERR "".$numberOfFragmentsB."\t(".sprintf("% .2f",100*$numberOfFragmentsB/$numberOfFragments)."%) "."\tbacterial fragments\n";
print STDERR "--------------------------------------------\n";


for(my $i=0;$i<=$#arrayofFilesbact;$i++){
  print STDERR "".$arrayofFilesbactToExtract[$i]."\t(".sprintf("% .2f",100*$arrayofFilesbactToExtract[$i]/$numberOfFragments)."%) "."\tbactaminant fragments from file: ".$arrayofFilesbact[$i]."\n";
}





########################
#                      #
#  Calling fragSim     #
#                      #
########################

#SELECTING ENDOGENOUS FRAGMENTS
if ($#arrayofFilesendo != -1 && $numberOfFragmentsE>0) {
  my $firstCmdFrag=1;
  if ($multiCellMode) {
  }else{
    $multiCellModeMAX=1;
  }

  for(my $i=0;$i<=$multiCellModeMAX;$i++) {
    if ($multiCellMode) {
      my @at_ = @{ $arrayOfEndoCellFasta[ $i ] };
      $numberOfFragmentsE1 = $arrayFragForEachCellE1[$i];
      $numberOfFragmentsE2 = $arrayFragForEachCellE2[$i];
      $arrayofFilesendo[0] =   $at_[0];
      if ($diploidMode) {
	$arrayofFilesendo[1] = $at_[1];
      }
    }

    if ($diploidMode) {
      my $cmd1="".$fragsim." -tag e1_".$i." -n ".$numberOfFragmentsE1;

      $cmd1 .= " -m ".$minsize." ";
      $cmd1 .= " -M ".$maxsize." ";

      if (defined $misince) {
	$cmd1 .= " --comp ".$misince." ";
	$cmd1 .= " --dist ".$distmis." ";
      }

      if (defined $methyl) {
	$cmd1 .= " --case ";
      }


      if (defined $filefragsize) {
	$cmd1 .= " -s ".$filefragsize." ";
      } else {
	if (defined $filefragfreqsize) {
	  $cmd1 .= " -f ".$filefragfreqsize." ";
	} else {
	  if (defined $loc) {
	    $cmd1 .= " --loc ".$loc." --scale ".$scale." ";
	  } else {
	    $cmd1 .= " -l ".$fraglength;
	  }
	}
      }


      if($uniq){
	$cmd1 .= " -uniq ";
      }

      $cmd1 .= "  ".$arrayofFilesendo[0]." |gzip ";
      if($firstCmdFrag){
	$cmd1 .= " > ";
	$firstCmdFrag=0;
      }else{
	$cmd1 .= " >> "
      }
      $cmd1 .= "  ".$outputprefix.".e.fa.gz";

      runcmd($cmd1);

      my $cmd2="".$fragsim." -tag e2_".$i." -n ".$numberOfFragmentsE2;

      $cmd2 .= " -m ".$minsize." ";
      $cmd2 .= " -M ".$maxsize." ";

      if (defined $methyl) {
	$cmd2 .= " --case ";
      }

      if (defined $filefragsize) {
	$cmd2 .= " -s ".$filefragsize." ";
      } else {
	if (defined $filefragfreqsize) {
	  $cmd2 .= " -f ".$filefragfreqsize." ";
	} else {
	  if (defined $loc) {
	    $cmd2 .= " --loc ".$loc." --scale ".$scale." ";
	  } else {
	    $cmd2 .= " -l ".$fraglength." ";
	  }
	}
      }

      if($uniq){
	$cmd2 .= " -uniq ";
      }

      $cmd2.=" ".$arrayofFilesendo[1]." |gzip >> ".$outputprefix.".e.fa.gz";
      runcmd($cmd2);

      #haploid mode
    } else {

      my $cmd1="".$fragsim." -tag e -n ".$numberOfFragmentsE1;

      $cmd1 .= " -m ".$minsize." ";
      $cmd1 .= " -M ".$maxsize." ";



      if (defined $methyl) {
	$cmd1 .= " --case ";
      }




      if (defined $misince) {
	$cmd1 .= " --comp ".$misince." ";
	$cmd1 .= " --dist ".$distmis." ";
      }

      if (defined $filefragsize) {
	$cmd1 .= " -s ".$filefragsize." ";
      } else {
	if (defined $filefragfreqsize) {
	  $cmd1 .= " -f ".$filefragfreqsize." ";
	} else {
	  if (defined $loc) {
	    $cmd1 .= " --loc ".$loc." --scale ".$scale." ";
	  } else {
	    $cmd1 .= " -l ".$fraglength." ";
	  }
	}
      }
      $cmd1 .= "  ".$arrayofFilesendo[0]." | gzip > ".$outputprefix.".e.fa.gz";
      runcmd($cmd1);

    }
  }				#end for each cell
} else {
  my $cmd1="touch ".$outputprefix.".e.fa.gz";
  runcmd($cmd1);
}

#SELECTING CONTAMINANT FRAGMENTS
if ($#arrayofFilescont != -1 && $numberOfFragmentsC>0) {

  for (my $i=0;$i<=$#arrayofFilescont;$i++) {
    my $cmd1="".$fragsim." -tag c".($i+1)." -n ".$arrayofFilescontToExtract[$i];

    $cmd1 .= " -m ".$minsize." ";
    $cmd1 .= " -M ".$maxsize." ";

    if (defined $methyl) {
      $cmd1 .= " --case ";
    }

    if (defined $misincc) {
      $cmd1 .= " --comp ".$misincc." ";
      $cmd1 .= " --dist ".$distmis." ";
    }

    if (defined $filefragsize) {
      $cmd1 .= " -s ".$filefragsize." ";
    } else {
      if (defined $filefragfreqsize) {
	$cmd1 .= " -f ".$filefragfreqsize." ";
      }else{
	if (defined $loc) {
	  $cmd1 .= " --loc ".$loc." --scale ".$scale." ";
	} else {
	  $cmd1 .= " -l ".$fraglength." ";
	}
      }
    }
    if($uniq){
      $cmd1 .= " -uniq ";
    }
    $cmd1 .= "  ".$arrayofFilescont[$i]." | gzip ";
    if($i==0){
      $cmd1 .= " >  ";
    }else{
      $cmd1 .= " >> ";
    }

    $cmd1.="".$outputprefix.".c.fa.gz";
    runcmd($cmd1);
  }
} else {
  my $cmd1="touch ".$outputprefix.".c.fa.gz";
  runcmd($cmd1);
}

#SELECTING BACTERIAL FRAGMENTS
if ($#arrayofFilesbact != -1 && $numberOfFragmentsB>0) {
  for (my $i=0;$i<=$#arrayofFilesbact;$i++) {
    my $cmd1="".$fragsim." -tag b".($i+1)." -n ".$arrayofFilesbactToExtract[$i];

    $cmd1 .= " -m ".$minsize." ";
    $cmd1 .= " -M ".$maxsize." ";

    if (defined $methyl) {
      $cmd1 .= " --case ";
    }

    if (defined $misincb) {
      $cmd1 .= " --comp ".$misincb." ";
    }

    if (defined $filefragsize) {
      $cmd1 .= " -s ".$filefragsize." ";
    } else {
      if (defined $filefragfreqsize) {
	$cmd1 .= " -f ".$filefragfreqsize." ";
      }else{
	if (defined $loc) {
	  $cmd1 .= " --loc ".$loc." --scale ".$scale." ";
	} else {
	  $cmd1 .= " -l ".$fraglength." ";
	}
      }
    }
    if($uniq){
      $cmd1 .= " -uniq ";
    }
    $cmd1 .= "  ".$arrayofFilesbact[$i]." | gzip ";
    if($i==0){
      $cmd1 .= " >  ";
    }else{
      $cmd1 .= " >> ";
    }
    $cmd1 .= " ".$outputprefix.".b.fa.gz";
    runcmd($cmd1);
  }
} else {
  my $cmd1="touch ".$outputprefix.".b.fa.gz";
  runcmd($cmd1);
}


########################
#                      #
#  Calling deamSim     #
#                      #
########################
#endogenous
if( (defined $matfilenonmeth)
    ||
    (defined $matfilemeth)
    ||
    (defined $matfilee)
    ||
    (defined $briggse)
    ||
    (@mapdamagee) ){

  my $cmde="".$deamsim." ";
  if( (defined $matfilee) ){
    $cmde .= " -matfile ".$matfilee." ";
  }

  if( (defined $matfilenonmeth) ){
    $cmde .= " -matfilenonmeth ".$matfilenonmeth." ";
  }

  if( (defined $matfilemeth) ){
    $cmde .= " -matfilemeth    ".$matfilemeth." ";
  }

  if( (defined $briggse) ){
    $cmde .= " -damage ".$briggse." ";
  }

  if( (@mapdamagee) ){
    $cmde .= " -mapdamage ".$mapdamagee[0]." ".$mapdamagee[1]." ";
  }

  $cmde.=" ".$outputprefix.".e.fa.gz  | gzip > ".$outputprefix."_d.fa.gz";
  runcmd($cmde);
}else{

  #no deamination on the endogenous
  my $cmde= "gzip -c -d  ".$outputprefix.".e.fa.gz  | gzip > ".$outputprefix."_d.fa.gz";
  runcmd($cmde);

}


#bacteria
if( (defined $matfilenonmeth)
    ||
    (defined $matfilemeth)
    ||
    (defined $matfileb)
    ||
    (defined $briggsb)
    ||
    (@mapdamageb) ){

  my $cmdb="".$deamsim." ";
  if( (defined $matfileb) ){
    $cmdb .= " -matfile ".$matfileb." ";
  }


  if( (defined $matfilenonmeth) ){
    $cmdb .= " -matfilenonmeth ".$matfilenonmeth." ";
  }

  if( (defined $matfilemeth) ){
    $cmdb .= " -matfilemeth    ".$matfilemeth." ";
  }

  if( (defined $briggsb) ){
    $cmdb .= " -damage ".$briggsb." ";
  }

  if( (@mapdamageb) ){
    $cmdb .= " -mapdamage ".$mapdamageb[0]." ".$mapdamageb[1]." ";
  }

  $cmdb.=" ".$outputprefix.".b.fa.gz  |gzip >> ".$outputprefix."_d.fa.gz";
  runcmd($cmdb);
}else{
  #no deamination on the bacteria
  my $cmdb= "gzip -c -d  ".$outputprefix.".b.fa.gz  | gzip >> ".$outputprefix."_d.fa.gz";
  runcmd($cmdb);
}










#human cont.
if( (defined $matfilenonmeth)
    ||
    (defined $matfilemeth)
    ||
    (defined $matfilec)
    ||
    (defined $briggsc)
    ||
    (@mapdamagec) ){

  my $cmdc="".$deamsim." ";
  if( (defined $matfilec) ){
    $cmdc .= " -matfile ".$matfilec." ";
  }

  if( (defined $matfilenonmeth) ){
    $cmdc .= " -matfilenonmeth ".$matfilenonmeth." ";
  }

  if( (defined $matfilemeth) ){
    $cmdc .= " -matfilemeth    ".$matfilemeth." ";
  }

  if( (defined $briggsc) ){
    $cmdc .= " -damage ".$briggsc." ";
  }

  if( (@mapdamagec) ){
    $cmdc .= " -mapdamage ".$mapdamagec[0]." ".$mapdamagec[1]." ";
  }

  $cmdc.=" ".$outputprefix.".c.fa.gz  |gzip >> ".$outputprefix."_d.fa.gz";
  runcmd($cmdc);
}else{

  #no deamination on the contaminant
  my $cmdec= "gzip -c -d  ".$outputprefix.".c.fa.gz  | gzip >> ".$outputprefix."_d.fa.gz";
  runcmd($cmdec);

}

########################
#                      #
#  Calling adptSim     #
#                      #
########################

my $cmdad =     "".$adptsim." ";
$cmdad   .= " -f ".$adapterF." ";
$cmdad   .= " -s ".$adapterR." ";
$cmdad   .= " -l ".$readlength." ";
if($se){
  $cmdad .= " -arts ".$outputprefix."_a.fa";
}else{
  $cmdad .= " -artp ".$outputprefix."_a.fa";
}
$cmdad .= "  ".$outputprefix."_d.fa.gz";
runcmd($cmdad);

########################
#                      #
#  Calling     art     #
#                      #
########################

my $cmdsq="".$artprog." -ss ".$ss." -amp -na ";
if($se){
  $cmdsq .= "    ";
}else{
  $cmdsq .= " -p ";
}
$cmdsq .= " -i ".$outputprefix."_a.fa ";
$cmdsq .= " -l ".$readlength." ";
$cmdsq .= " -c 1 ";
$cmdsq .= " -qs  ".$qs." ";
$cmdsq .= " -qs2 ".$qs2." ";

$cmdsq .= " -o ".$outputprefix."_s";
runcmd($cmdsq);

my $cmdzip;

$cmdzip = "gzip -f ".$outputprefix."_a.fa";
runcmd($cmdzip);

if($se){
  $cmdzip = "gzip -f ".$outputprefix."_s.fq";
  runcmd($cmdzip);
  print STDERR "Single-end reads available here: ".$outputprefix."_s.fq.gz\n";
}else{
  $cmdzip = "gzip -f ".$outputprefix."_s1.fq";
  runcmd($cmdzip);
  $cmdzip = "gzip -f ".$outputprefix."_s2.fq";
  runcmd($cmdzip);

  print STDERR "Paired-end forward reads available here: ".$outputprefix."_s1.fq.gz\n";
  print STDERR "Paired-end reverse reads available here: ".$outputprefix."_s2.fq.gz\n";

}

print STDERR "Program finished successfully\n";
