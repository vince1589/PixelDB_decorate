#!/usr/bin/perl
use List::Util qw(shuffle);
use strict;

# Script that take a list of command, a split, an id an run them and a working directory

my $pwd = `pwd`;
chomp($pwd);

my $saut = -1;
my $id = -1;
my $list = "NA";
my $holo = "";
my $wdir = "NA";
my $random = 0;
my $test = 0;
foreach my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq "-sp") {$saut = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-id") {$id = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-list") {$list = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-folder") {$wdir = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-rand") {$random = 1;}
	if ($ARGV[$i] eq "-test") {$test = 1;}
}

unless(-e $list) {"Please input a list: -list -> $list\n";die;}
if ($wdir eq "NA") {"Please input a working dir: -folder\n";die;}

unless ($wdir =~ m/^\//) {$wdir = "$pwd/$wdir";}
#Because id can only start at 1, will add -1
$id = $id - 1;
print "This is ID=$id\n";

# Load list
open IN, "$list" or die "Can't open list -> $list\n";
my @list = <IN>;
close IN;
my $scalar = scalar(@list);
if ($random == 1) {
    print "I shuffle list with seed = 12\n";
    srand(12);
    @list = shuffle @list;

}
print "I have $scalar command\n";
my @split = ();
if ($saut == -1) {
	#print "I don't aggregate the list\n";
	foreach my $i (0..@list-1) {
		my @temp = ($list[$i]);
		push(@split,\@temp);
	}
} else {
	
	my $mod = int($scalar / $saut+0.5);
	print "I split the list in $saut parts, wich is composed of $mod\n";
	my $c = 0;
	my $push = 0;
	foreach my $i (0..@list-1) {
		push(@{$split[$c]},$list[$i]);
		++$push;
		if ($i != 0 && $push >= $mod) {
			++$c;
			
			$mod = int(($scalar-$i) / ($saut-$c)+0.5);
			#printf "I pushed:$push -> $i $c NM:$mod = ($scalar-$i) / ($saut-$c) = %d / %d = %f\n",($scalar-$i),($saut-$c),($scalar-$i) / ($saut-$c);
			$push = 0;

		}
	}
	#print "Last is $push\n";
}
#die;
# Run command

foreach my $i (0..@split-1) {
	#print "$i $id\n";
	if ($id != -1 && $id != $i) {next;}
	#print "HERE\n";
	#print "$split[$i]\n";
	foreach my $f (@{$split[$i]}) {
		chomp($f);
		chdir $wdir;

        # Run command
        print "$f\n";
	if ($test == 0) {system "$f";}
		
		
	}
}

print "Done!\n";










