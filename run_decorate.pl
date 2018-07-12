#!/usr/bin/perl;

use List::Util qw(shuffle);

use strict;



#/home/grigoryanlab/home/cmack2357/scripts/structgen/programs/decorate --targetPDB /home/ironfs/scratch/grigoryanlab/vfrap/decorate/test/apo_3fdl.pdb --outDir /home/ironfs/scratch/grigoryanlab/vfrap/decorate/test/ --maxRMSD 1.1 --persLen 15 --scoreMethod 1 --scoreCut -20 --maxNumSegs 8 --maxNumFreeSegs 6 --first 0 --last 20 --interCutBB 4.0 --maxRMSDBB 0.7 --interCutS 0.03 --maxRMSDS 0.8 --minRatio 0.1
my $pdb = "NA";
my $pepch = "NA";
my $wrkdir = "NA";
my $min = 0;
my $max = 10;;
my $overwrite = 0;
my $random = 0;
my $dry = 0;


foreach my $i (0..@ARGV-1) {
  if ($ARGV[$i] eq "-pdb") {$pdb = $ARGV[$i+1];} #PDB to run with all the peptide chain
  if ($ARGV[$i] eq "-pepch") {$pepch = $ARGV[$i+1];} #Specify the peptide chain (comma separated)
  if ($ARGV[$i] eq "-wrkdir") {$wrkdir = $ARGV[$i+1];} # Where you want to work
  if ($ARGV[$i] eq "-min") {$min = $ARGV[$i+1];} #First TERMs to use (in the ranked order by abundance)
  if ($ARGV[$i] eq "-max") {$max = $ARGV[$i+1];} #Last TERMs to use (in the ranked order by abundance)
  if ($ARGV[$i] eq "-random") {$random = 1;} #Do you want to shuffle ranked order (usefull if you have multiple job that are running same subset)
  
  #Overwrite level (0 == If job is deemed to be running will not run it, 1 == Will overwrite job that are running, 2 == job that are done)
  #Sometime job will crash and system will think it is still running, so use overwritting == 1 when re-lauching job
  if ($ARGV[$i] eq "-overwrite") {$overwrite = $ARGV[$i+1];} 
  
  if ($ARGV[$i] eq "-dry") {$dry = 1;print "Dry = 1\n";} #Will not run any decorate
  
}

print "Usage\n-pdb\t\tPDB file to decorate\n-pepch\t\tChain of the peptide\n-wrkdir\t\tWorking directory\n-min\t\tmin TERMs\n-max\t\tmax Terms\n\n\n";

if ($pdb eq "NA") {die;}


#Gen PDB name
my $pdbname = "";
if ($pdb =~ m/([A-Za-z0-9_]+)\.pdb$/) {$pdbname = $1;} else {die;}

print "PDBNAME=$pdbname\n";

# If no wrkdir, use default
if ($wrkdir eq "NA") {
  my $path = "/home/ironfs/scratch/grigoryanlab/vfrap/decorate/";
  $wrkdir = $path . $pdbname ."/";
  print "No wrkdir provided, will use default\n";
}

#Make wrkdir
print "WRKDIR= $wrkdir\n";
unless (-d "$wrkdir") {mkdir $wrkdir;print "Creating wrkdir\n";}
unless (-d "$wrkdir/run") {mkdir "$wrkdir/run";print "Creating run\n";} #Where the TERMs completion are created
unless (-d "$wrkdir/info") {mkdir "$wrkdir/info";print "Creating info\n";} #Info is where job status are stored (ran/completed)
unless (-d "$wrkdir/results") {mkdir "$wrkdir/results";print "Creating results\n";} #This where the results are stored (overlap with crystal structure peptide), more on that

#Find peptide chain if from PixelDB

if ($pepch eq "NA") {print "Please specify peptide chain (comma separated), will exit\n";die;}
my @pepch_sp = split(",",$pepch);
#Copy Holo
system "cp $pdb $wrkdir/";

#Generate apo form
my $apopdb = "$wrkdir/$pdbname\_apo\.pdb";
if ((! -e $apopdb) | ($overwrite == 1)) {
  print "Creating apo file:$wrkdir/$pdbname\_apo\.pdb\n";
  open OUT, ">$wrkdir/$pdbname\_apo\.pdb" or die;
  open IN, "$pdb";
  while(my $l = <IN>) {
    my $tonext = 0;
    foreach my $ch (@pepch_sp) {
      if ($l =~ m/^ATOM................ $pepch /) {$tonext += 1;} #Next if peptide chain (Print
    }
    if ($tonext != 0) {next;}
    print OUT "$l";
  }
  close IN;
  close OUT;
}

if ($dry == 1) {exit();}

my @arr = ($min..$max);
if ($random == 1) {my @arr = shuffle @arr;}

#Load all done
my %alldone = load_all_done("$wrkdir/info/all_done.dat"); #I don't remember how this file is created lol, but it keep track of what is already run (faster than looking up for file in a folder)

#Run each TERMs example (one at the time... this might not be the most effecient, check with Craig maybe)
foreach my $i (@arr) {
  if (exists $alldone{$i}) {next;}
  if ($overwrite < 1) {if (-e "$wrkdir/info/running_$i.dat") {next;}}
  if ($overwrite < 2) {if (-e "$wrkdir/info/done_$i.dat") {next;}}
  my $upperbound = $i + 1;
  
  #Create file that say that I'm running something
  system "touch $wrkdir/info/running_$i.dat";
  my $out = `/home/grigoryanlab/home/cmack2357/scripts/structgen/programs/old/decorate --targetPDB $apopdb --outDir $wrkdir/run/ --maxRMSD 1.1 --persLen 15 --scoreMethod 0 --scoreCut -20 --maxNumSegs 8 --maxNumFreeSegs 6 --first $i --last $upperbound --interCutBB 4.0 --maxRMSDBB 0.7 --interCutS 0.03 --maxRMSDS 0.8 --minRatio 0.1`;
  print "$out\n";
  
  #If it is done, delete file that says it is running and create file that say it is done
  if ($out =~ m/PROCEDURE COMPLETE/) {
    system "rm $wrkdir/info/running_$i.dat";
    system "touch $wrkdir/info/done_$i.dat";
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
















