#!/usr/bin/perl

my $max = 50000;
my $step = 10000;

#Short run
#$max = 100;
#$step = 50;


open IN, "ecrs_list.dat" or die; #List of file to run
my $path = "/home/ifs-users/vfrap/";

open LST, ">job/run_decorate_list.dat" or die;
open DEC, ">job/run_analyse_list.dat" or die;
open CON, ">job/run_concat_list.dat" or die;

while (my $l = <IN>) {
  chomp($l);
  #Get clusters
  my $clust = "NA";
  my $ch = "NA";
  if ($l =~ m/_(.)_(\d+)_\d+\.pdb/) {$clust = $2;$ch = $1;}
  
  #Create file with dry run
  #This will create holo, folder,apo file... so that process don't try to write at the same time
  `perl $path/decorate/script/run_decorate.pl -pdb $path/PixelDB/clusters/$clust/$l -min 0 -max 0 -pepch $ch -dry`;
   
  for my $i (0..$max) {
    my $nowmin = $i*$step;
    my $nowmax = ($i+1)*$step-1;
    if ($nowmin > $max) {last;}
    if ($nowmax > $max) {$nowmax = $max;}
    print LST "perl $path/decorate/script/run_decorate.pl -pdb $path/PixelDB/clusters/$clust/$l -min $nowmin -max $nowmax -pepch $ch\n";
  }
  my $name = $l;
  $name =~ s/\.pdb//;
  print DEC "python $path/decorate/script/analyse_decorate_folder.py $path/PixelDB/clusters/$clust/$l $ch $path/decorate/results/$name.pickle /home/ironfs/scratch/grigoryanlab/vfrap/decorate/$name/run/*.pdb\n";
  print CON "perl $path/decorate/script/concat_everything.pl /home/ironfs/scratch/grigoryanlab/vfrap/decorate/$name/\n";
  
  
  
  #last;
  #perl script/run_decorate.pl -pdb ../PixelDB/clusters/99/1JD5_A_B_99_1.pdb -max 10

}
close IN;
close LST;
close DEC;
close CON;
