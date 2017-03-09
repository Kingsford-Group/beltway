#!/usr/bin/perl
use strict;

while(<>){
	chomp $_;
	if($_ =~ /.*cycloLen([0-9]*)_.*User time \(seconds\):\s+([0-9\.]*)/){
		print "$1\t$2\n";
	}
	else{
		print "-1\t-1\n";
	}
}
