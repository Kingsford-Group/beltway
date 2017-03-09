#!/usr/bin/perl
use strict;

while(<>){
	chomp $_;
	if($_ =~ /.*cycloLen([0-9]*)_.*or m:ss\):\s+([0-9]*):([0-9]*):([0-9]*)/){
		print "$1\t".($4+($3 + ($2 * 60) * 60))."\n";
	}
	elsif($_ =~ /.*cycloLen([0-9]*)_.*or m:ss\):\s+([0-9]*):([0-9\.]*)/){
		print "$1\t".($3 + ($2 * 60))."\n";
	}
}
