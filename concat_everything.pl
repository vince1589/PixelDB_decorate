#!/usr/bin/perl

use strict;

#Script that concat a bunch of file so it easier to read


foreach my $fold (@ARGV) {

  #Chdir to folder (everything will be relative)
  chdir "$fold";
  
  my %alldone = load_all_done("./info/all_done.dat");
  foreach my $f (glob "./info/done_*.dat") {
    if ($f =~ m/done_(\d+)\.dat/) {$alldone{$1} += 1;} else {die;}
  }
  
  #Write output
  my @arr = sort {$a <=> $b} keys %alldone;
  print "@arr\n";
  my $str = "";
  my $buf = $arr[0];
  foreach my $i (1..@arr) {
    if (($arr[$i] - $arr[$i-1] != 1) | ($i == scalar(@arr))) {
      if ($arr[$i-1] eq $buf) {} else {
      $buf .= "-$arr[$i-1]";
      }
      $str .= "$buf,";
      $buf = "$arr[$i]";
    }
  }
  $str =~ s/,$//;
  open OUT, ">./info/all_done.dat" or die;
  print OUT "$str\n";
  close OUT;
  
  foreach my $f (glob "./info/done_*.dat") {
    system "rm $f";
  }
  
}



sub load_all_done() {
   my $f = $_[0];
  #Concat all the done
  my %alldone;
  # Load predone
  open IN, "$f";
  while (my $l = <IN>){
    chomp($l);
    my @sp = split(/,/,$l);
    foreach my $s (@sp) {
      if ($s =~ m/^(\d+)-(\d+)/) {
        my $one = $1;my $two = $2;
        for my $i ($one..$two) {$alldone{$i} +=1;}
      } else {
        $alldone{$s} += 1;
      }
    }
  }
  return(%alldone);
}
