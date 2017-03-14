#!/usr/bin/perl
use strict;

foreach my $file(`ls data/float_simul*/*.wPeptide.txt`){
	my %list;
	chomp $file;
	$file =~ /.*cycloLen([0-9]*)_.*/;
	my $len = $1;

	foreach my $line(`cat $file`){
		chomp $line;
		my @spl = split(/\s+/,$line);
		my @AA = split(/\-/,$spl[1]);
		if($list{$spl[0]} == 1 && $len != scalar(@AA)){
			print $file."\t".$spl[0]."\t$len\t".scalar(@AA)."\n";
			last;
		}
		$list{$spl[0]} = 1;
	}
}
