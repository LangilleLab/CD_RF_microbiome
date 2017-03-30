use strict;
use warnings;

if ( scalar @ARGV != 3 )	{	

	print "This script takes in 3 arguments:\n";
	print "1) A VCF containing the SNPs (with rsIDs) or interest\n";
	print "2) A file containing the sample order output by shapeit (which when I ran it was the same for each chromosome)\n";
	print "3) The path to the impute2 output files, specifically the \"haps\" files that are named in this format: RSID_impute2_output_haps (e.g. rs16967103_impute2_output_haps)\n";
	print "A new VCF containing the imputed genotypes for each SNP will be returned\n";

	die "usage: perl impute2vcf.pl (VCF of rsIDs of interest) (file containing shapeit output order) (impute2 output folder) > imputed.vcf\n";	
}

my $known_vcf = $ARGV[0];
my $sample_order = $ARGV[1];
my $impute2_output = $ARGV[2];

my %known = (); # hash containing rsIDs as keys and the corresponding full line of the VCF as values (the first few columns are used in the output VCF)

open( 'KNOWN' , '<' , $known_vcf ) or die "cant open KNOWN $known_vcf\n";
while( <KNOWN> )	{

	# Skip header lines
	my $first = substr( $_ , 0 , 1 );
	if ( $first eq "#" )	{	next	}

	# Split by tabs and newlines since VCFs are tab-delimited
	my @s = split( '[\t\n]', $_ );

	# Add key/value pair of rsID and full VCF line in array
	$known{$s[2]} = \@s;
}

my @impute_order = (); # array containing the sample orders
my $lc = 0;
open( 'SAMPLES' , '<' , $sample_order ) or die "cant open SAMPLES $sample_order\n";
while( <SAMPLES> )	{

	# Skip first two lines
	if ( $lc < 2 )	{	++$lc;	next	}

	my @s = split( '\s+' , $_ );

	# add sample name to array
	push( @impute_order , $s[0] );

} close( 'SAMPLES' );

# Print out header
print "##fileformat=VCFv4.2\n";
my @header = ("#CHROM" , qw(POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT) , @impute_order); 
my $h = join( "	" , @header);
print "$h\n";

# Open impute2 output directory
opendir(DIR, $impute2_output ) or die "cannot open DIR $impute2_output\n";

# Read in all files that end in "impute2_output_haps"
my @f = grep(/_impute2_output_haps$/ ,readdir(DIR));

# Close directory
closedir(DIR);

# Loop through all these files
foreach my $f ( @f )    {

	# Split names by "_" and figure out rsID
    my @file_info = split( '_' , $f );
    my $rsID = $file_info[0];

    # Figure out full path to file
	my $path = $impute2_output . "/" . $f;

	# Open haps file and loop through line by line
	open( 'IMPUTED' , '<', $path ) or die "cant open IMPUTED $path\n";
	while( <IMPUTED> )	{

		my @s = split( '\s+' , $_ );

		# If row is not for focal rsID then skip line
		my @info = split( ':' , $s[1] );
		if ( $rsID ne $info[0] )	{	next	}	

		# Otherwise continue and get original VCF info
		# Print out first 8 columns
		my @original = @{$known{$info[0]}};
		my @start = ( splice( @original, 0, 8 ), "GT" );
		my $start = join( "	" , @start );
		print "$start";

		# Read in all imputed genotype calls
		my @calls = splice( @s , 5 );
		
		# Counter for index of calls array
		my $i = 0;

		# Loop through all individuals in order
		foreach my $ind ( @impute_order )	{

			# Figure out base1 and base2 for individual and print out genotype in VCF format
			my $b1 = $calls[$i];

			# Add 1 to get the next genotype
			++$i;
			my $b2 = $calls[$i];

			# Print genotype
			my $code = $b1 ."/". $b2;
			print "	$code";

			# Add 1 so next loop is for next individual
			++$i; 
			
		}
		# Print newline after all genotypes are outputted
		print "\n";

	} close( 'IMPUTED' );
}
