#!/usr/bin/perl
use strict;

my $given_fname = shift;
my $solution_fname = shift;

if($given_fname eq "" || $solution_fname eq ""){
	die("Usage ./check.pl <test_file> <reference_file>\n");
}

our @given;
open GFILE, $given_fname or die("$given_fname: $!\n");
while(<GFILE>){
	chomp $_;
	my @spl = split(/\-/,$_);
	foreach my $aa(@spl){
		push @given, $aa;
	}
}

our @solution;
open SFILE, $solution_fname or die("$solution_fname: $!\n");
while(<SFILE>){
        chomp $_;
        my @spl = split(/\-/,$_);
        foreach my $aa(@spl){
                push @solution, $aa;
        }
}

sub check_position_direction{
	my $start_in_given = shift;
	my $direction = shift;
	my $j = $start_in_given;
	for(my $i=0;$i<scalar(@solution);$i++){
                return 0 if $solution[$i] ne $given[$j];
		if($direction eq "reverse"){
			$j = ($j==0)?(scalar(@given)-1):($j-1);
		}else{
			$j = ($j + 1) % scalar(@given);
		}
	}
	return 1;
}

for my $i(0...scalar(@given)-1){
	if($solution[0] eq $given[$i]){
		if(check_position_direction($i) == 1){
			print "The solution maps to the answer at position $i\n";
                        exit(0);
		}elsif(check_position_direction($i,"reverse") == 1){
			print "The solution maps to the answer reversed at position $i\n";
			exit(0);
		}
	}	
}

print "No mapping found\n";


