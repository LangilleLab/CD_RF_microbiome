#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;

if ( scalar @ARGV != 2 )	{	
	print "This script takes in 2 arguments, a VCF and an output directory. It will output 1 VCF for each chromsome (i.e. for each unique value in first column of VCF)\n";
	die "usage: perl splitVCFbyChr.pl (vcf) (out dir)\n";	
}

my $vcf = $ARGV[0];
my $outdir = $ARGV[1];

system("mkdir -p $outdir");

my %fh = (); # chromosomes are keys and and filehandles are values
my @h = (); # array containing all headerlines

open( 'VCF' , '<' , $vcf ) or die "cant open VCF $vcf\n";
while (<VCF> )	{
	
	my $first = substr( $_ , 0 , 1 );
	if ( $first eq "#" )	{	 
		# first character is # keep in "@h" variable and skip
		push( @h , $_ );
		next;	
	}		

	my @s = split( '[\t\n]' , $_ );

	my $chr = $s[0];

	if ( ! exists $fh{$chr} )	{

		# first line of chr so make new output file
		my $out = basename( $vcf );
		$out =~ s/\.vcf$/.$chr.vcf/;
		$out = $outdir . "/" . $out;

		# open up output filehandle
		open( $fh{$chr} , '>' , $out ) or die "cant create $out\n";

		# print out full header to vcf
		my $header = join("" , @h );
		print {$fh{$chr}} "$header";
	}	

	# print line to VCF for that particular chr
	print {$fh{$chr}} "$_";

} close( 'VCF' );
