#!/usr/bin/perl


# This script is a user-unfriendly way to create the diploid genomes files from a methylation map (see my $fileMethyl= "/[path to file]/chrom20.WGBS.methylation.filtered.bed";)
# it will open the first argument which is the genome file. It will set certain bases as methylated and others as not methylated
# 

use strict;
use warnings;
use Data::Dumper;

sub selectRandom{
  my ($t,$num)=@_;
  my @array;

  for (my $i=0;$i<$num;$i++) {
    $array[$i]=0;
  }

  my $succ=0;
  while ($succ<$t) {
    my $i=int(rand($num));
    if ($array[$i]==0) {
      $array[$i]=1;
      $succ++;
    }
  }
  return @array;
}

#my @t=selectRandom(20,100);

#print Dumper(@t);
#die;

my %methylCoord;
my %methylCoordFound;

my $fileMethyl= "/[path to file]/chrom20.WGBS.methylation.filtered.bed";
#my $fileMethyl= "chrom20.WGBS.methylation.filtered.bed";
open(FILEM, $fileMethyl) || die "cannot open pipe to $fileMethyl";

while(my $line = <FILEM>){
  chomp($line);
  #print $line."\n";
  my @fields=split("\t",$line);
  #if( ($fields[10] >=50) &&
  #    ($fields[5] eq "+")    )      {
    #print $line."\n";
    #push(@methylCoord,$fields[1]);
  if($fields[10] > 0){
    $methylCoord{$fields[1]}     =$fields[10];
    $methylCoordFound{$fields[1]}=0;
  }

  #}
}
close(FILEM);
#die;
warn "done";
#warn "Will label ".$#methylCoord." sites as methylated\n";
#die;

my @arrayFD;

for(my $i=0;$i<100;$i++){
  open($arrayFD[ $i ] , ">chr20_".$i.".fa");
}

my $file= $ARGV[0];
open(FILE, "gunzip -c $file |") || die "cannot open pipe to $file";
my $coord=0;
my $lastChr="";
my $line = <FILE>;
#print $line;

for(my $i=0;$i<100;$i++){
  print { $arrayFD[ $i ] } $line;
}

#for(my $i=0;$i<100;$i++){
#  close($arrayFD[ $i ] );
#}

my $lastWasCpgM=0;
my @arrayWhichMethyl;
my $lineSeen=0;
while($line = <FILE>){
  chomp($line);
  my $l=length($line);
  #warn "saw line $line\n";
  $lineSeen++;
  for(my $i=0;$i<$l;$i++){

    my $currentChr = substr($line,$i,1);
    #print $coord."\t".$currentChr."\t".$lastChr."\t".$lastWasCpgM."\n";

    if($lastChr     eq "C" &&
       $currentChr  eq "G"){

      if( (exists $methylCoord{$coord-1}) ){

	#print $coord."\t".$lastChr."\t".$currentChr."\t1\n";
	$methylCoordFound{$coord-1}=1;


	@arrayWhichMethyl = selectRandom($methylCoord{$coord-1},100);
	#die Dumper(@arrayWhichMethyl);
	for(my $j=0;$j<100;$j++){
	  if($arrayWhichMethyl[ $j ] == 1){
	    print { $arrayFD[ $j ] } lc($lastChr);
	  }else{
	    print { $arrayFD[ $j ] } uc($lastChr);
	  }
	}

	#print lc($lastChr);
	$lastChr = $currentChr;
	$lastWasCpgM=1;
      }else{

	#if( (exists $methylCoord{$coord}) ){
	#  $methylCoordFound{$coord}=1;
	#  print lc($lastChr);
	#  $lastChr = lc($currentChr);
	#  die $coord;
	#}else{
	#print $lastChr;
	for(my $j=0;$j<100;$j++){
	  print { $arrayFD[ $j ] } $lastChr;
	}

	$lastChr =    $currentChr;
	#}
      }

    }else{
      #print $lastChr;
      if( $lastWasCpgM ){
	for(my $j=0;$j<100;$j++){
	  if($arrayWhichMethyl[ $j ] == 1){
	    print { $arrayFD[ $j ] } lc($lastChr);
	  }else{
	    print { $arrayFD[ $j ] } uc($lastChr);
	  }
	}
	$lastWasCpgM=0;
      }else{

	for(my $j=0;$j<100;$j++){
	  print { $arrayFD[ $j ] } $lastChr;
	}

      }

      $lastChr = $currentChr;

    }


    $coord++;

    if($i==0 && $lineSeen!=1){
       for(my $i=0;$i<100;$i++){
	 print { $arrayFD[ $i ] } "\n";
       }
     }

  }#end for each char


  #print "\n";
}
close(FILE);

for(my $j=0;$j<100;$j++){
  print { $arrayFD[ $j ] } $lastChr."\n";
}

for(my $i=0;$i<100;$i++){
  close($arrayFD[ $i ] );
}


foreach my $c (keys %methylCoordFound){
  if($methylCoordFound{$c} == 0){
    warn "no ".$c."\n";
  }else{
    warn "yes ".$c."\n";
  }
}
