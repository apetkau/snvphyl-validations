#!/usr/bin/env perl
# Purpose
# Given two positions tables of variants, determines the number of unique and non-unique positions

use strict;
use warnings;

use Set::Scalar;

# read positions file a Set::Scalar object
sub read_pos
{
	my ($file) = @_;

	my %results;

	open(my $fh,"<$file") or die "Could not open $file\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	$results{'header'} = $header;
	$results{'positions'} = Set::Scalar->new;
	while(my $line = readline($fh))
	{
		chomp($line);
		my @tokens = split(/\t/,$line);

		my ($chrom,$position,$status,$reference) = @tokens;
		die "Error with line $line, position invalid\n" if (not defined $position or $position !~ /\d+/);

		# ignore anything not valid
		next if ($status ne 'valid');

		$results{'positions'}->insert($line);
	}
	close($fh);

	return \%results;
}

my $usage = "$0 [variants1.tsv] [variants2.tsv]\n";

die "error:\n$usage" if (@ARGV != 2);

my ($variants1_f,$variants2_f) = @ARGV;

my $variants1 = read_pos($variants1_f);
my $variants2 = read_pos($variants2_f);

# must have same genomes and same order of genomes
if ($variants1->{'header'} ne $variants2->{'header'})
{
	die "Error: headers did not match\n";
}

my $variants1_positions = $variants1->{'positions'};
my $variants2_positions = $variants2->{'positions'};

# set operations
my $intersection = $variants1_positions * $variants2_positions;
my $uniq_variants1 = $variants1_positions - $variants2_positions;
my $uniq_variants2 = $variants2_positions - $variants1_positions;

print "$variants1_f\t$variants2_f\tIntersection\tUnique-$variants1_f\tUnique-$variants2_f\n";
print $variants1_positions->size,"\t",$variants2_positions->size,"\t",$intersection->size,
	"\t",$uniq_variants1->size,"\t",$uniq_variants2->size,"\n";
