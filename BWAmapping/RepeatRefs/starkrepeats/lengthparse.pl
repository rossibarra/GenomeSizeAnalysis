#!/usr/bin/perl
use strict; use warnings;

my $input = $ARGV[0];
my $output = $ARGV[1];

open(BLAST, "<$input");
open(FILT, ">$output");

while(<BLAST>){
	chomp;
	my ($query, $sub, $perc, $alin) = split("\t");
	if($alin > 30){
		if($perc > 80){
			print FILT "$_\n";
		}
	}
}
