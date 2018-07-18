#!/usr/bin/perl

use strict;

my $bin_path = "/media/vince/Postdoc/side_project/PixelDB_decorate/";

my $folder = $ARGV[0]; #Full path!!!

#Make full path if not
my $pwd = `pwd`;
chomp($pwd);
unless ($folder =~ m/^\//) {$folder = "$pwd/$folder";}


#Get in folder of results
chdir $folder;

#Create results folder
mkdir "pickle_result_files";

#Make the binding mode
unless (-d "$folder/binding_pose") {
  print "python $bin_path/merge_frag_fast.py $folder\n";
  system "python $bin_path/merge_frag_fast.py $folder";
}

#Make list of frament

unless (-e "./Fragment.lst") {
  system "ls -1 ./*nres* ./binding_pose/*nres* > Fragment.lst";
}
#Get PDB to run
my @glob = glob "./bound_complex/*.pdb";

foreach my $f (@glob) {
  print "$f\n";
  #Find chain
  my $ch = "";
  my $name;
  if ($f =~ m/bound_complex\/(.*)\.pdb/) {
    $name = $1;
    my @sp = split("_",$1);
    $ch = $sp[2];
  } else {
    die;
  }
  unless (-e "./pickle_result_files/$name.pk") {
    my $cmd = "python $bin_path/analyse_decorate_folder.py $f $ch ./pickle_result_files/$name.pk Fragment.lst";
    print "$cmd\n";
    system "$cmd";
    #last;
  }
  
  #python ../analyse_decorate_folder.py ~/PostDoc/PixelDB/gitPixelDB/PixelDB/clusters/1/1AGC_A_C_1_2.pdb C 1A1M_decorate_results/pickle_file/1AGC_A_C.pk FragMent.lst
}

